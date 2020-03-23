#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "fpga_codec.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"
#ifdef ENABLE_FPGA
#include <utils/lcd.h>
#include "dma_common.h"
#endif
#include <pthread.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <poll.h>
#include <sys/time.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define	MEM_16G		(1ULL << 34)
#define BATCH_SIZE  1000
#define BATCH_LINE_LIMIT	16384
#define TIMEOUT     BATCH_SIZE*100*1000      // Nanoseconds
#define MIN(x,y)    ((x < y)? x : y)
typedef fpga_pci_data_t fpga_pci_conn;
#define NUM_FPGA_THREADS	1
/* Theory on probability and scoring *ungapped* alignment
 *
 * s'(a,b) = log[P(b|a)/P(b)] = log[4P(b|a)], assuming uniform base distribution
 * s'(a,a) = log(4), s'(a,b) = log(4e/3), where e is the error rate
 *
 * Scale s'(a,b) to s(a,a) s.t. s(a,a)=x. Then s(a,b) = x*s'(a,b)/log(4), or conversely: s'(a,b)=s(a,b)*log(4)/x
 *
 * If the matching score is x and mismatch penalty is -y, we can compute error rate e:
 *   e = .75 * exp[-log(4) * y/x]
 *
 * log P(seq) = \sum_i log P(b_i|a_i) = \sum_i {s'(a,b) - log(4)}
 *   = \sum_i { s(a,b)*log(4)/x - log(4) } = log(4) * (S/x - l)
 *
 * where S=\sum_i s(a,b) is the alignment score. Converting to the phred scale:
 *   Q(seq) = -10/log(10) * log P(seq) = 10*log(4)/log(10) * (l - S/x) = 6.02 * (l - S/x)
 *
 *
 * Gap open (zero gap): q' = log[P(gap-open)], r' = log[P(gap-ext)] (see Durbin et al. (1998) Section 4.1)
 * Then q = x*log[P(gap-open)]/log(4), r = x*log[P(gap-ext)]/log(4)
 *
 * When there are gaps, l should be the length of alignment matches (i.e. the M operator in CIGAR)
 */


#ifdef __cplusplus
extern "C" {
#endif

#define QUEUESIZE 1000
typedef struct{
	bseq1_t **seqs;
	mem_alnreg_v ** regs;
	mem_chain_v ** chains;
	int64_t num;
	int64_t starting_read_id;

	// For loading
	int last_entry;
} queue_t;


typedef struct {
	queue_t **buf;
	long head, tail;
	int full, empty;
	pthread_mutex_t *mut;
	pthread_cond_t *notFull, *notEmpty;
} queue;









typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	const mem_pestat_t *pes;
	smem_aux_t **aux;
	bseq1_t *seqs;
	int64_t n_processed;
	queue *queue1;
} worker_t;





typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	const mem_pestat_t *pes;
	queue *queue1;
} worker2_t;

typedef struct {
	worker_t *w_master;     // Location of master with all sequences
	int tid;                // Thread id

	// Sequences processed by any thread will be all seqs starting from tid*BATCH_SIZE;
	// next batch to be processed will be opt->n_threads*BATCH_SIZE

} worker_slave_t;


typedef struct {
	queue *q1;      // Queue for stage 1 - 2 ( worker1_MT  |   q1   | fpga_worker)
	queue *q2;      // Queue for stage 2 - 3 ( fpga_worker |   q2   | worker2_MT)
	worker_t * w;
	pthread_mutex_t *seedex_mut;
	int tid;
} queue_coll;       // Collection of queues










uint64_t total_seeds = 0;

queue *queueInit (void)
{
	queue *q;

	q = (queue *)malloc (sizeof (queue));
	if (q == NULL) return (NULL);

	q->empty = 1;
	q->full = 0;
	q->head = 0;
	q->tail = 0;
	q->mut = (pthread_mutex_t *) malloc (sizeof (pthread_mutex_t));
	pthread_mutex_init (q->mut, NULL);
	q->notFull = (pthread_cond_t *) malloc (sizeof (pthread_cond_t));
	pthread_cond_init (q->notFull, NULL);
	q->notEmpty = (pthread_cond_t *) malloc (sizeof (pthread_cond_t));
	pthread_cond_init (q->notEmpty, NULL);

	q->buf = (queue_t**) malloc(QUEUESIZE * sizeof(queue_t*));
	
	return (q);
}
void queueDelete (queue *q)
{
	pthread_mutex_destroy (q->mut);
	free (q->mut);      
	pthread_cond_destroy (q->notFull);
	free (q->notFull);
	pthread_cond_destroy (q->notEmpty);
	free (q->notEmpty);
	free (q->buf);
	free (q);
}
void queueAdd (queue *q, queue_t* in)
{
	q->buf[q->tail] = in;
	q->tail++;
	if (q->tail == QUEUESIZE)
		q->tail = 0;
	if (q->tail == q->head)
		q->full = 1;
	q->empty = 0;

	return;
}
void queueDel (queue *q, queue_t **out)
{
	*out = q->buf[q->head];

	q->head++;
	if (q->head == QUEUESIZE)
		q->head = 0;
	if (q->head == q->tail)
		q->empty = 1;
	q->full = 0;

	return;
}

int fpga_verbose = 0;

	// FPGA Write Buffers


	const uint32_t write_buffer_read_entry_length = 64;    // Each read entry is 64 bytes
	const uint32_t write_buffer_seed_entry_length = 32;    // Each seed entry is 32 bytes
	const uint32_t write_buffer_num_entries = 256;
	const uint32_t write_buffer_capacity = 256 * 64;
	const uint32_t max_write_buffer_index = 256 * 64;


	uint64_t fpga_mem_write_offset = 0;

	uint32_t vled;
	uint32_t vdip;
	
	pci_bar_handle_t bw_pci_bar_handle;

	fpga_pci_conn * fpga_pci_local;


	struct timeval s2_waitq1_st, s2_waitq1_et;
	int total_s1_waitq1_time = 0;





static const bntseq_t *global_bns = 0; // for debugging only

mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1; o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;

	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1<<30;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

/***************************
 * Collection SA invervals *
 ***************************/

#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)


static smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
	return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{	
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

static void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a)
{
	int i, k, x = 0, old_n;
	int start_width = 1;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	a->mem.n = 0;
	// first pass: find all SMEMs
	while (x < len) {
		if (seq[x] < 4) {
			x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i) {
				bwtintv_t *p = &a->mem1.a[i];
				int slen = (uint32_t)p->info - (p->info>>32); // seed length
				if (slen >= opt->min_seed_len)
					kv_push(bwtintv_t, a->mem, *p);
			}
		} else ++x;
	}
	// second pass: find MEMs inside a long SMEM
	old_n = a->mem.n;
	for (k = 0; k < old_n; ++k) {
		bwtintv_t *p = &a->mem.a[k];
		int start = p->info>>32, end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
		for (i = 0; i < a->mem1.n; ++i)
			if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len)
				kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
	}
	// third pass: LAST-like
	if (opt->max_mem_intv > 0) {
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				if (1) {
					bwtintv_t m;
					x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
					if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
				} else { // for now, we never come to this block which is slower
					x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
					for (i = 0; i < a->mem1.n; ++i)
						kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
				}
			} else ++x;
		}
	}
	// sort
	ks_introsort(mem_intv, a->mem.n, a->mem.a);
}


#include "kbtree.h"

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

// return 1 if the seed is merged into the chain
static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const mem_seed_t *p, int seed_rid)
{
	int64_t qend, rend, x, y;
	const mem_seed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (seed_rid != c->rid) return 0; // different chr; request a new chain
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) && p->rbeg >= l_pac) return 0; // don't chain if on different strand
	x = p->qbeg - last->qbeg; // always non-negtive
	y = p->rbeg - last->rbeg;
	if (y >= 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
		if (c->n == c->m) {
			c->m <<= 1;
			c->seeds = realloc(c->seeds, c->m * sizeof(mem_seed_t));
		}
		c->seeds[c->n++] = *p;
		return 1;
	}
	return 0; // request to add a new chain
}

int mem_chain_weight(const mem_chain_t *c)
{
	int64_t end;
	int j, w = 0, tmp;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->qbeg >= end) w += s->len;
		else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
		end = end > s->qbeg + s->len? end : s->qbeg + s->len;
	}
	tmp = w; w = 0;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->rbeg >= end) w += s->len;
		else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
		end = end > s->rbeg + s->len? end : s->rbeg + s->len;
	}
	w = w < tmp? w : tmp;
	return w < 1<<30? w : (1<<30)-1;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn)
{
	int i, j;
	for (i = 0; i < chn->n; ++i) {
		mem_chain_t *p = &chn->a[i];
		err_printf("* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
		for (j = 0; j < p->n; ++j) {
			bwtint_t pos;
			int is_rev;
			pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= p->seeds[j].len - 1;
			err_printf("\t%d;%d;%d,%ld(%s:%c%ld)", p->seeds[j].score, p->seeds[j].len, p->seeds[j].qbeg, (long)p->seeds[j].rbeg, bns->anns[p->rid].name, "+-"[is_rev], (long)(pos - bns->anns[p->rid].offset) + 1);
		}
		err_putchar('\n');
	}
}

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf)
{
	int i, b, e, l_rep;
	int64_t l_pac = bns->l_pac;
	mem_chain_v chain;
	kbtree_t(chn) *tree;
	smem_aux_t *aux;

	kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);

	aux = buf? (smem_aux_t*)buf : smem_aux_init();
	mem_collect_intv(opt, bwt, len, seq, aux);
	for (i = 0, b = e = l_rep = 0; i < aux->mem.n; ++i) { // compute frac_rep
		bwtintv_t *p = &aux->mem.a[i];
		int sb = (p->info>>32), se = (uint32_t)p->info;
		if (p->x[2] <= opt->max_occ) continue;
		if (sb > e) l_rep += e - b, b = sb, e = se;
		else e = e > se? e : se;
	}
	l_rep += e - b;
	for (i = 0; i < aux->mem.n; ++i) {
		bwtintv_t *p = &aux->mem.a[i];
		int step, count, slen = (uint32_t)p->info - (p->info>>32); // seed length
		int64_t k;
		// if (slen < opt->min_seed_len) continue; // ignore if too short or too repetitive
		step = p->x[2] > opt->max_occ? p->x[2] / opt->max_occ : 1;
		for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
			mem_chain_t tmp, *lower, *upper;
			mem_seed_t s;
			int rid, to_add = 0;
			s.rbeg = tmp.pos = bwt_sa(bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
			s.qbeg = p->info>>32;
			s.score= s.len = slen;
			rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
			if (rid < 0) continue; // bridging multiple reference sequences or the forward-reverse boundary; TODO: split the seed; don't discard it!!!
			if (kb_size(tree)) {
				kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
				if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid)) to_add = 1;
			} else to_add = 1;
			if (to_add) { // add the seed as a new chain
				tmp.n = 1; tmp.m = 4;
				tmp.seeds = calloc(tmp.m, sizeof(mem_seed_t));
				tmp.seeds[0] = s;
				tmp.rid = rid;
				tmp.is_alt = !!bns->anns[rid].is_alt;
				kb_putp(chn, tree, &tmp);
			}
		}
	}
	if (buf == 0) smem_aux_destroy(aux);

	kv_resize(mem_chain_t, chain, kb_size(tree));

	#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
	#undef traverse_func

	for (i = 0; i < chain.n; ++i) chain.a[i].frac_rep = (float)l_rep / len;
	if (bwa_verbose >= 4) printf("* fraction of repetitive seeds: %.3f\n", (float)l_rep / len);

	kb_destroy(chn, tree);
	return chain;
}

/********************
 * Filtering chains *
 ********************/

#define chn_beg(ch) ((ch).seeds->qbeg)
#define chn_end(ch) ((ch).seeds[(ch).n-1].qbeg + (ch).seeds[(ch).n-1].len)

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, mem_chain_t, flt_lt)

int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a)
{
	int i, k;
	kvec_t(int) chains = {0,0,0}; // this keeps int indices of the non-overlapping chains
	if (n_chn == 0) return 0; // no need to filter
	// compute the weight of each chain and drop chains with small weight
	for (i = k = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		c->first = -1; c->kept = 0;
		c->w = mem_chain_weight(c);
		if (c->w < opt->min_chain_weight) free(c->seeds);
		else a[k++] = *c;
	}
	n_chn = k;
	ks_introsort(mem_flt, n_chn, a);
	// pairwise chain comparisons
	a[0].kept = 3;
	kv_push(int, chains, 0);
	for (i = 1; i < n_chn; ++i) {
		int large_ovlp = 0;
		for (k = 0; k < chains.n; ++k) {
			int j = chains.a[k];
			int b_max = chn_beg(a[j]) > chn_beg(a[i])? chn_beg(a[j]) : chn_beg(a[i]);
			int e_min = chn_end(a[j]) < chn_end(a[i])? chn_end(a[j]) : chn_end(a[i]);
			if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
				int li = chn_end(a[i]) - chn_beg(a[i]);
				int lj = chn_end(a[j]) - chn_beg(a[j]);
				int min_l = li < lj? li : lj;
				if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
					large_ovlp = 1;
					if (a[j].first < 0) a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
					if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
						break;
				}
			}
		}
		if (k == chains.n) {
			kv_push(int, chains, i);
			a[i].kept = large_ovlp? 2 : 3;
		}
	}
	for (i = 0; i < chains.n; ++i) {
		mem_chain_t *c = &a[chains.a[i]];
		if (c->first >= 0) a[c->first].kept = 1;
	}
	free(chains.a);
	for (i = k = 0; i < n_chn; ++i) { // don't extend more than opt->max_chain_extend .kept=1/2 chains
		if (a[i].kept == 0 || a[i].kept == 3) continue;
		if (++k >= opt->max_chain_extend) break;
	}
	for (; i < n_chn; ++i)
		if (a[i].kept < 3) a[i].kept = 0;
	for (i = k = 0; i < n_chn; ++i) { // free discarded chains
		mem_chain_t *c = &a[i];
		if (c->kept == 0) free(c->seeds);
		else a[k++] = a[i];
	}
	return k;
}

/******************************
 * De-overlap single-end hits *
 ******************************/

#define alnreg_slt2(a, b) ((a).re < (b).re)
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

#define alnreg_slt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define alnreg_hlt(a, b)  ((a).score > (b).score || ((a).score == (b).score && ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

#define alnreg_hlt2(a, b) ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash2, mem_alnreg_t, alnreg_hlt2)

#define PATCH_MAX_R_BW 0.05f
#define PATCH_MIN_SC_RATIO 0.90f

int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b, int *_w)
{
	int w, score, q_s, r_s;
	double r;
	if (bns == 0 || pac == 0 || query == 0) return 0;
	assert(a->rid == b->rid && a->rb <= b->rb);
	if (a->rb < bns->l_pac && b->rb >= bns->l_pac) return 0; // on different strands
	if (a->qb >= b->qb || a->qe >= b->qe || a->re >= b->re) return 0; // not colinear
	w = (a->re - b->rb) - (a->qe - b->qb); // required bandwidth
	w = w > 0? w : -w; // l = abs(l)
	r = (double)(a->re - b->rb) / (b->re - a->rb) - (double)(a->qe - b->qb) / (b->qe - a->qb); // relative bandwidth
	r = r > 0.? r : -r; // r = fabs(r)
	if (bwa_verbose >= 4)
		printf("* potential hit merge between [%d,%d)<=>[%ld,%ld) and [%d,%d)<=>[%ld,%ld), @ %s; w=%d, r=%.4g\n",
			   a->qb, a->qe, (long)a->rb, (long)a->re, b->qb, b->qe, (long)b->rb, (long)b->re, bns->anns[a->rid].name, w, r);
	if (a->re < b->rb || a->qe < b->qb) { // no overlap on query or on ref
		if (w > opt->w<<1 || r >= PATCH_MAX_R_BW) return 0; // the bandwidth or the relative bandwidth is too large
	} else if (w > opt->w<<2 || r >= PATCH_MAX_R_BW*2) return 0; // more permissive if overlapping on both ref and query
	// global alignment
	w += a->w + b->w;
	w = w < opt->w<<2? w : opt->w<<2;
	if (bwa_verbose >= 4) printf("* test potential hit merge with global alignment; w=%d\n", w);
	bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, pac, b->qe - a->qb, query + a->qb, a->rb, b->re, &score, 0, 0);
	q_s = (int)((double)(b->qe - a->qb) / ((b->qe - b->qb) + (a->qe - a->qb)) * (b->score + a->score) + .499); // predicted score from query
	r_s = (int)((double)(b->re - a->rb) / ((b->re - b->rb) + (a->re - a->rb)) * (b->score + a->score) + .499); // predicted score from ref
	if (bwa_verbose >= 4) printf("* score=%d;(%d,%d)\n", score, q_s, r_s);
	if ((double)score / (q_s > r_s? q_s : r_s) < PATCH_MIN_SC_RATIO) return 0;
	*_w = w;
	return score;
}

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a)
{
	int m, i, j;
	if (n <= 1) return n;
	ks_introsort(mem_ars2, n, a); // sort by the END position, not START!
	for (i = 0; i < n; ++i) a[i].n_comp = 1;
	for (i = 1; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		if (p->rid != a[i-1].rid || p->rb >= a[i-1].re + opt->max_chain_gap) continue; // then no need to go into the loop below
		for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
			mem_alnreg_t *q = &a[j];
			int64_t or_, oq, mr, mq;
			int score, w;
			if (q->qe == q->qb) continue; // a[j] has been excluded
			or_ = q->re - p->rb; // overlap length on the reference
			oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (or_ > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
				if (p->score < q->score) {
					p->qe = p->qb;
					break;
				} else q->qe = q->qb;
			} else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w)) > 0) { // then merge q into p
				p->n_comp += q->n_comp + 1;
				p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
				p->sub = p->sub > q->sub? p->sub : q->sub;
				p->csub = p->csub > q->csub? p->csub : q->csub;
				p->qb = q->qb, p->rb = q->rb;
				p->truesc = p->score = score;
				p->w = w;
				q->qb = q->qe;
			}
		}
	}
	for (i = 0, m = 0; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	n = m;
	ks_introsort(mem_ars, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i-1].score && a[i].rb == a[i-1].rb && a[i].qb == a[i-1].qb)
			a[i].qe = a[i].qb;
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	return m;
}

typedef kvec_t(int) int_v;

static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z)
{ // similar to the loop in mem_chain_flt()
	int i, k, tmp;
	tmp = opt->a + opt->b;
	tmp = opt->o_del + opt->e_del > tmp? opt->o_del + opt->e_del : tmp;
	tmp = opt->o_ins + opt->e_ins > tmp? opt->o_ins + opt->e_ins : tmp;
	z->n = 0;
	kv_push(int, *z, 0);
	for (i = 1; i < n; ++i) {
		for (k = 0; k < z->n; ++k) {
			int j = z->a[k];
			int b_max = a[j].qb > a[i].qb? a[j].qb : a[i].qb;
			int e_min = a[j].qe < a[i].qe? a[j].qe : a[i].qe;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].sub == 0) a[j].sub = a[i].score;
					if (a[j].score - a[i].score <= tmp && (a[j].is_alt || !a[i].is_alt))
						++a[j].sub_n;
					break;
				}
			}
		}
		if (k == z->n) kv_push(int, *z, i);
		else a[i].secondary = z->a[k];
	}
}

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id)
{
	int i, n_pri;
	int_v z = {0,0,0};
	if (n == 0) return 0;
	for (i = n_pri = 0; i < n; ++i) {
		a[i].sub = a[i].alt_sc = 0, a[i].secondary = a[i].secondary_all = -1, a[i].hash = hash_64(id+i);
		if (!a[i].is_alt) ++n_pri;
	}
	ks_introsort(mem_ars_hash, n, a);
	mem_mark_primary_se_core(opt, n, a, &z);
	for (i = 0; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		p->secondary_all = i; // keep the rank in the first round
		if (!p->is_alt && p->secondary >= 0 && a[p->secondary].is_alt)
			p->alt_sc = a[p->secondary].score;
	}
	if (n_pri >= 0 && n_pri < n) {
		kv_resize(int, z, n);
		if (n_pri > 0) ks_introsort(mem_ars_hash2, n, a);
		for (i = 0; i < n; ++i) z.a[a[i].secondary_all] = i;
		for (i = 0; i < n; ++i) {
			if (a[i].secondary >= 0) {
				a[i].secondary_all = z.a[a[i].secondary];
				if (a[i].is_alt) a[i].secondary = INT_MAX;
			} else a[i].secondary_all = -1;
		}
		if (n_pri > 0) { // mark primary for hits to the primary assembly only
			for (i = 0; i < n_pri; ++i) a[i].sub = 0, a[i].secondary = -1;
			mem_mark_primary_se_core(opt, n_pri, a, &z);
		}
	} else {
		for (i = 0; i < n; ++i)
			a[i].secondary_all = a[i].secondary;
	}
	free(z.a);
	return n_pri;
}

/*********************************
 * Test if a seed is good enough *
 *********************************/

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200

#define MEM_HSP_COEF 1.1f
#define MEM_MINSC_COEF 5.5f
#define MEM_SEEDSW_COEF 0.05f

int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_seed_t *s)
{
	int qb, qe, rid;
	int64_t rb, re, mid, l_pac = bns->l_pac;
	uint8_t *rseq = 0;
	kswr_t x;

	if (s->len >= MEM_SHORT_LEN) return -1; // the seed is longer than the max-extend; no need to do SW
	qb = s->qbeg, qe = s->qbeg + s->len;
	rb = s->rbeg, re = s->rbeg + s->len;
	mid = (rb + re) >> 1;
	qb -= MEM_SHORT_EXT; qb = qb > 0? qb : 0;
	qe += MEM_SHORT_EXT; qe = qe < l_query? qe : l_query;
	rb -= MEM_SHORT_EXT; rb = rb > 0? rb : 0;
	re += MEM_SHORT_EXT; re = re < l_pac<<1? re : l_pac<<1;
	if (rb < l_pac && l_pac < re) {
		if (mid < l_pac) re = l_pac;
		else rb = l_pac;
	}
	if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN) return -1; // the seed seems good enough; no need to do SW

	rseq = bns_fetch_seq(bns, pac, &rb, mid, &re, &rid);
	x = ksw_align2(qe - qb, (uint8_t*)query + qb, re - rb, rseq, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, 0);
	free(rseq);
	return x.score;
}

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn, mem_chain_t *a)
{
	double min_l = opt->min_chain_weight? MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
	int i, j, k, min_HSP_score = (int)(opt->a * min_l + .499);
	if (min_l > MEM_SEEDSW_COEF * l_query) return; // don't run the following for short reads
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		for (j = k = 0; j < c->n; ++j) {
			mem_seed_t *s = &c->seeds[j];
			s->score = mem_seed_sw(opt, bns, pac, l_query, query, s);
			if (s->score < 0 || s->score >= min_HSP_score) {
				s->score = s->score < 0? s->len * opt->a : s->score;
				c->seeds[k++] = *s;
			}
		}
		c->n = k;
	}
}

/****************************************
 * Construct the alignment from a chain *
 ****************************************/

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

#define MAX_BAND_TRY  2




void encode_read_data(char * read_seq, int length, uint32_t read_id, uint8_t *ar) {
	//memset(ar, 0, ar_size*sizeof(uint8_t));
	
	int i = 0; 

	if(bwa_verbose >= 18){
		printf("[READ]0,%u,%d,",read_id,length);
		int j = 0;
		for(j=0;j<length;j++){
			if(read_seq[j] == 0) {
				printf("A");
			}   
			else if(read_seq[j] == 1){
				printf("C");
			} 
			else if(read_seq[j] == 2){
				printf("G");
			} 
			else if(read_seq[j] == 3){
				printf("T");
			}
			else {
				printf("N");
			}
		}
		printf(",\n");

	}

	ar[0] = 1;
	memcpy(ar+1,&read_id,4);
	ar[5] = (uint8_t)length; 

	// Use a counter for the number of bits already encoded. Never go beyond 8.
	// Encode the char and increment the counter by 3. 
	// If the number is less than 3, (3-x) is the bits to be put in the 
	// previous entry and x bits in the new entry

	int count = 0;
	int buffer_index = 6;
	uint8_t read_char = 0;

	for(i=0;i<length;i++){
		if(read_seq[i] == 0) {
			read_char = 0;
		}   
		else if(read_seq[i] == 1){
			read_char = 1;
		} 
		else if(read_seq[i] == 2){
			read_char = 2;
		} 
		else if(read_seq[i] == 3){
			read_char = 3;
		}
		else{ 
			read_char = 4;
		}
		if(count == 0){
			count += 3;
			ar[buffer_index] = read_char;
		}
		else{
			count += 3;
			count %= 8;
			if(count < 3){
				ar[buffer_index] = ar[buffer_index] | (read_char<<(8-(3-count)));
				buffer_index++;
				ar[buffer_index] = ar[buffer_index] | (read_char>>(3-count));
			}
			else{
				ar[buffer_index] = ar[buffer_index] | (read_char<<(count-3));
			}       
			
		}
	}

}


	

	int encode_seed_data(uint64_t abs_chain_ref_beg, uint64_t abs_seed_ref_beg, uint64_t abs_chain_ref_end, uint32_t is_reverse, uint32_t seed_score,uint32_t seed_read_beg, uint32_t seed_read_end, uint32_t rel_seed_ref_beg, uint32_t rel_seed_ref_end, uint8_t *ar, int ar_index) {
	 
		int ar_offset = 1;
		//memset(ar, 0, ar_size*sizeof(uint32_t));

		if(bwa_verbose >= 18){
			printf("[SEED]1,%lu,%lu,%lu,%u,%u,%u,%u,%u,%u,",abs_chain_ref_beg,abs_seed_ref_beg,abs_chain_ref_end,is_reverse,seed_score,seed_read_beg,seed_read_end,rel_seed_ref_beg,rel_seed_ref_end);

			if(seed_read_beg != 0) {
				printf("1,1,");
			}
			else {
				printf("0,0,");
			}
			if(seed_read_end != 101) {
				printf("1,0,\n");
			}
			else {
				printf("0,0,\n");
			}
			fflush(stdout);
		}

		if(ar_index == 0){
			ar[0] = 3;
		}
		else{
			ar[0] = 1;
		}

		memcpy(&ar[ar_offset],&abs_chain_ref_beg,8);
		ar_offset += 8;

		memcpy(&ar[ar_offset],&abs_seed_ref_beg,8);
		ar_offset += 8;

		memcpy(&ar[ar_offset],&abs_chain_ref_end,8);
		ar_offset += 8;

		ar[ar_offset] = is_reverse;
		ar_offset++;

		ar[ar_offset] = (uint8_t)(seed_score & 0xFF);
		ar_offset++;
		
		ar[ar_offset] = (uint8_t)(seed_read_beg & 0xFF);
		ar_offset++;

		ar[ar_offset] = (uint8_t)(seed_read_end & 0xFF);
		ar_offset++;

		ar[ar_offset] = (uint8_t)(rel_seed_ref_beg & 0xFF);
		ar_offset++;

		ar[ar_offset] = (uint8_t)(rel_seed_ref_end & 0xFF);
		ar_offset++;


		uint8_t ext_int= 0;


		if(seed_read_end != 101) {
			ext_int = 1;
			//push_to_lsb(ar+3,ar_size-3,0x0,1);
			//push_to_lsb(ar+3,ar_size-3,0x1,1);
		}
		else {
			//push_to_lsb(ar+3,ar_size-3,0x0,1);
			//push_to_lsb(ar+3,ar_size-3,0x0,1);
		}
		ext_int = (ext_int << 2);

		if(seed_read_beg != 0) {
			ext_int = (ext_int | 3);
			//push_to_lsb(ar+3,ar_size-3,0x1,1);
			//push_to_lsb(ar+3,ar_size-3,0x1,1);
		}
		else {
			//push_to_lsb(ar+3,ar_size-3,0x0,1);
			//push_to_lsb(ar+3,ar_size-3,0x0,1);
		}

		ar[ar_offset] = ext_int;
		ar_offset++;
		
		if(ar_index == 0) {
			return 256;
		} 
		else {
			return 0;
		}

	}









void fetch_rmaxs(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av, int64_t* rmax0, int64_t* rmax1){
	int i; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], max = 0;


	if (c->n == 0) return;
		
			 
	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len;
	}
	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;


	if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
		if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else rmax[0] = l_pac;
	}
		if (bwa_verbose >= 10) {
				printf("*** FPGA : rmax[0] :%ld, rmax[1]: %ld   \n",rmax[0],rmax[1]);
		}       

		*rmax0 = rmax[0];
		*rmax1 = rmax[1];
}



/*
	Pack 8 byte encoded string (size of len) to 3 byte string (8 chars)
*/
void f_8to3(const char *a, int len, uint8_t *b)
{
	int i, j;
	uint32_t *p = (uint32_t*)b;

	for (i = 0; i < len; i += 8)
	{
		int offset = 0;
		*p = 0;
		for (j = 0; j < 8; j++)
		{
			*p |= (uint32_t) a[i+j] << offset;
			offset += 3;
		}
		//printf("%d: %lo\n", i, *p);
		p = (uint32_t *)((char*)p + 3);

	}

}

void f_3to8(const char *a, int len, char *b)
{
	int i, j;
	const uint32_t *p = (const uint32_t*)a;
	
	for (i = 0; i < len; i += 8)
	{
		int offset = 0;
		for (j = 0; j < 8; j++)
		{
			b[i+j] = (*p & (0x7 << offset)) >> (offset);
			offset += 3;
		}
		p = (const uint32_t *)((const char*)p + 3);
	}
}

struct SeedExPackageGen
{
	SeedExPackageGen(): has_next(false) {}

	void * new_input (LineParams params, char * query, char * target, union SeedExLine* buf)
	{
		char payload_buf[72];
		buf->ty1.preamble = PACKET_START;
		buf->ty1.params = params;
		query_ptr = query;
		target_ptr = target;
		padding = params.tlen - params.qlen;
		qlen = params.qlen;
		tlen = params.tlen;

		assert(padding >= 0);

		// Query
		// memset(buf->ty1.payload1, 0, 27);
		if (params.qlen + padding < 72) /* |query------|padding---|nul---| */
		{
			memcpy(payload_buf, query, params.qlen);
			memset(payload_buf + params.qlen, C_PADDING, padding);
			memset(payload_buf + params.qlen + padding, C_NULL, 72 - params.qlen - padding);
			query_ptr += params.qlen;
			padding = 0;
			has_next = false;
		}
		else if (params.qlen < 72) /* |query------|padding---| */
		{
			memcpy(payload_buf, query, params.qlen);
			memset(payload_buf + params.qlen, C_PADDING, 72 - params.qlen);
			query_ptr += params.qlen;
			padding -= 72 - params.qlen;
			has_next = true;
		}
		else /* |query-------------------| */
		{
			memcpy(payload_buf, query, 72);
			query_ptr += 72;
			has_next = true;
		}

		// printf("Query:\t"); for (int i = 0; i < 72; ++i) printf("%d", payload_buf[i]); putchar('\n');
		f_8to3(payload_buf, 72, buf->ty1.payload1);
		qlen -= query_ptr - query;

		// Target
		// memset(buf->ty1.payload2, 0, 27);
		if (params.tlen < 72) /* |target------|nul---| */
		{
			memcpy(payload_buf, target, params.tlen);
			memset(payload_buf + params.tlen, C_NULL, 72 - params.tlen);
			target_ptr += params.tlen;
		}
		else /* |target-------------------|| */
		{
			memcpy(payload_buf, target, 72);
			target_ptr += 72;
		}
		// printf("Target:\t"); for (int i = 0; i < 72; ++i) printf("%d", payload_buf[i]); putchar('\n');
		f_8to3(payload_buf, 72, buf->ty1.payload2);
		tlen -= target_ptr - target;

		if (!has_next) {
			buf->ty1.preamble = PACKET_END;
		}
	}

	void * next(union SeedExLine* buf)
	{
		if (!has_next) return nullptr;

		char payload_buf[88];
		char * orig_q = query_ptr, * orig_t = target_ptr;
		buf->ty0.preamble = PACKET_MIDDLE;
		// memset(buf->ty0.payload, 0xffff, 63);

		// Query
		if (qlen + padding < 84) /* |query------|padding---|nul---| */
		{
			memcpy(payload_buf, query_ptr, qlen);
			memset(payload_buf + qlen, C_PADDING, padding);
			memset(payload_buf + qlen + padding, C_NULL, 84 - qlen - padding);
			query_ptr += qlen;
			padding = 0;
			has_next = false;
		}
		else if (qlen < 84) /* |query------|padding---| */
		{
			memcpy(payload_buf, query_ptr, qlen);
			memset(payload_buf + qlen, C_PADDING, 84 - qlen);
			query_ptr += qlen;
			padding -= 84 - qlen;
			assert(padding >= 0);
			has_next = true;
		}
		else /* |query-------------------| */
		{
			memcpy(payload_buf, query_ptr, 84);
			query_ptr += 84;
			has_next = true;
		}
		// printf("Query(o):\t"); for (int i = 0; i < 84; ++i) printf("%d", payload_buf[i]); putchar('\n');
		// printf("Query(w):\t    "); for (int i = 0; i < 80; ++i) printf("%d", payload_buf[i]); putchar('\n');
		f_8to3(payload_buf, 80, buf->ty0.payload);
		qlen -= query_ptr - orig_q;

		// Target
		memcpy(payload_buf, payload_buf + 80, 4);
		char * buf_target_ptr = &payload_buf[4];
		if (tlen < 84) /* |target------|nul---| */
		{
			memcpy(buf_target_ptr, target_ptr, tlen);
			memset(buf_target_ptr + tlen, C_NULL, 84 - tlen);
			target_ptr += tlen;
		}
		else /* |target-------------------|| */
		{
			memcpy(buf_target_ptr, target_ptr, 84);
			target_ptr += 84;
		}
		// printf("Target(w):\t"); for (int i = 0; i < 88; ++i) printf("%d", payload_buf[i]); putchar('\n');
		f_8to3(payload_buf, 88, &buf->ty0.payload[30]);
		tlen -= target_ptr - orig_t;

		if (!has_next) {
			buf->ty1.preamble = PACKET_END;
		}
	}

	char * query_ptr;
	char * target_ptr;
	int padding;
	bool has_next;
	int qlen, tlen;
};


void mem_chain2aln_to_fpga(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av, int64_t rmax0, int64_t rmax1, fpga_data_out_t* fpga_result, std::vector<union SeedExLine>& write_buffer1, std::vector<union SeedExLine*>& write_buffer_entry1, std::vector<union SeedExLine>& write_buffer2, std::vector<union SeedExLine*>& write_buffer_entry2, std::vector<struct extension_meta_t>& extension_meta)
{
	int i, k, rid, aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	uint64_t *srt;
	SeedExPackageGen gen;

	if (c->n == 0) return;
		// FPGA : Write read data into write_buffer
		if (bwa_verbose >= 10) {
				int j;
				printf("*** FPGA : Seeing Read Query:   "); for (j = 0; j < l_query; ++j) putchar("ACGTN"[(int)query[j]]); putchar('\n');
		}

		rmax[0] = rmax0;
		rmax[1] = rmax1;

	// retrieve the reference sequence
	rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
	// if(bwa_verbose >= 15){
	//     rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
	// }
	//
	//int64_t k3 = 0;
	//assert(c->rid == rid);

	// srt = malloc(c->n * 8);
	// for (i = 0; i < c->n; ++i)
	// 	srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
	// ks_introsort_64(c->n, srt);

	// for (k = c->n - 1; k >= 0; --k) {
	for (k = 0; k < c->n; ++k) {
		mem_alnreg_t *a;
		// s = &c->seeds[(uint32_t)srt[k]];
		s = &c->seeds[k];

		if(s->qbeg == 0 && ((s->qbeg + s->len) == l_query)){
			a = kv_pushp(mem_alnreg_t, *av);
			memset(a, 0, sizeof(mem_alnreg_t));
			a->w = aw[0] = aw[1] = opt->w;
			a->score = a->truesc = -1;
			a->rid = c->rid;
			a->score = a->truesc = s->len * opt->a;
			a->qb = 0;
			a->rb = s->rbeg;
			a->qe = l_query;
			a->re = s->rbeg + s->len;
			a->seedlen0 = s->len;
			a->frac_rep = c->frac_rep;
		}
		else{
			// uint32_t is_reverse = (s->rbeg >= (l_pac)) ? 1 : 0; 

			int seq_id = write_buffer_entry1.size();
			a = kv_pushp(mem_alnreg_t, *av);
			memset(a, 0, sizeof(mem_alnreg_t));
			a->w = aw[0] = aw[1] = opt->w;
			a->score = a->truesc = -1;
			a->rid = c->rid;

			if(bwa_verbose >= 15){
				int j = 0;
				printf("[REFERENCE] %ld,",rmax[1] - rmax[0]); for (j = 0; j < (rmax[1] - rmax[0]) ; ++j) putchar("ACGTN"[(int)rseq[j]]); putchar('\n');
			}
			//*ar_index = encode_seed_data(rmax[0], s->rbeg, rmax[1], is_reverse,s->len,s->qbeg, (s->qbeg + s->len), (uint32_t)(s->rbeg - rmax[0]), (uint32_t)(s->rbeg - rmax[0] + s->len), *write_buffer + *write_buffer_index, *ar_index);

			if (s->qbeg) { // left extension
				uint8_t *rs, *qs;
				int qle, tle, gtle, gscore;
				qs = malloc(s->qbeg);
				for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
				tmp = s->rbeg - rmax[0];
				rs = malloc(tmp);
				for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];

				write_buffer1.push_back({});
				write_buffer_entry1.push_back(&write_buffer1.back());
				gen.new_input({seq_id, s->qbeg, tmp, s->len * opt->a, aw[0]}, (char*)qs, (char*)rs, &write_buffer1.back());
				while (gen.has_next)
				{
					write_buffer1.push_back({});
					gen.next(&write_buffer1.back());
				}
				// ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
				// if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
				fpga_result->fpga_entry_present = 1;
				free(qs); free(rs);
			} else {
				a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;
				write_buffer_entry1.push_back(nullptr);
			}

			if (s->qbeg + s->len != l_query) { // right extension
				int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
				qe = s->qbeg + s->len;
				re = s->rbeg + s->len - rmax[0];
				assert(re >= 0);

				write_buffer2.push_back({});
				write_buffer_entry2.push_back(&write_buffer2.back());
				gen.new_input({seq_id, l_query - qe, rmax[1] - rmax[0] - re, sc0, aw[1]}, (char*)query + qe, (char*)rseq + re, &write_buffer2.back());
				while (gen.has_next)
				{
					write_buffer2.push_back({});
					gen.next(&write_buffer2.back());
				}
				//ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
				// if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
				fpga_result->fpga_entry_present = 1;
			} else {
				a->qe = l_query, a->re = s->rbeg + s->len;
				write_buffer_entry2.push_back(nullptr);
			}

			extension_meta.back().seed_id = k;
			extension_meta.push_back(extension_meta.back());

			a->seedlen0 = s->len;
			a->frac_rep = c->frac_rep;

			if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);

			total_seeds++;
		}

	}
	free(rseq);
}


void postprocess_alnreg (const mem_opt_t *opt, int l_query, const mem_chain_t *c, mem_alnreg_v *av0, mem_alnreg_v *av)
{
	int i, k;

	uint64_t *srt;
	srt = malloc(c->n * 8);
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
	ks_introsort_64(c->n, srt);

	for (k = c->n - 1; k >= 0; --k) {
		mem_alnreg_t *a;
		auto s = &c->seeds[(uint32_t)srt[k]];                    // Select seed with best score first within a chain

		for (i = 0; i < av->n; ++i) { // test whether extension has been made before
			mem_alnreg_t *p = &av->a[i];
			int64_t rd;
			int qd, w, max_gap;
			if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
			if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
			// qd: distance ahead of the seed on query; rd: on reference
			qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
			max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
			w = max_gap < p->w? max_gap : p->w; // bounded by the band width
			if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
			// similar to the previous four lines, but this time we look at the region behind
			qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
			max_gap = cal_max_gap(opt, qd < rd? qd : rd);
			w = max_gap < p->w? max_gap : p->w;
			if (qd - rd < w && rd - qd < w) break;
		}
		if(bwa_verbose >= 18){
				printf("(FPGA) i = %d,k = %d, av_size = %zu\n",i,k,av->n);
		}
		if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
			if (bwa_verbose >= 4)
				printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
					   k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
			for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
				const mem_seed_t *t;
				if (srt[i] == 0) continue;
				t = &c->seeds[(uint32_t)srt[i]];
				if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
				if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
				if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
			}
			if (i == c->n) { // no overlapping seeds; then skip extension
				srt[k] = 0; // mark that seed extension has not been performed
				continue;
			}
			if (bwa_verbose >= 4)
				printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
		}

		a = kv_pushp(mem_alnreg_t, *av);
		memcpy(a, &av0->a[(uint32_t)srt[k]], sizeof(mem_alnreg_t));
	}

	free(srt);
}


void fpga_func_model(const mem_opt_t *opt, std::vector<union SeedExLine>& load_buf, std::vector<union SeedExLine*>& idx, std::vector<union SeedExLine>& results)
{
	int i,j;
	char buf[168];
	for (auto p : idx)
	{
		if (!p) continue;
		int qlen = p->ty1.params.qlen;
		int tlen = p->ty1.params.tlen;
		int w = p->ty1.params.w;
		int init_score = p->ty1.params.init_score;
		int seq_id = p->ty1.params.seq_id;
		char *query = calloc(qlen, sizeof(char));
		char *target = calloc(tlen, sizeof(char));
		char *query_ptr = query;
		char *target_ptr = target;
		assert(p->ty1.preamble == PACKET_START || p->ty1.preamble == PACKET_END);

		// decode first line
		f_3to8(p->ty1.payload1, 72, buf);
		memcpy(query, buf, MIN(72, qlen));
		query_ptr += MIN(72, qlen);

		f_3to8(p->ty1.payload2, 72, buf);
		memcpy(target, buf, MIN(72, tlen));
		target_ptr += MIN(72, tlen);

		if (p->ty1.preamble == PACKET_START) {
			// while (target_ptr - target < tlen + 1)
			for (union SeedExLine *line = (union SeedExLine *)p + 1; ; ++line)
			{
				assert(line->ty0.preamble == PACKET_MIDDLE || line->ty0.preamble == PACKET_END);
				f_3to8(line->ty0.payload, 168, buf);
				if (qlen > query_ptr - query) {
					memcpy(query_ptr, buf, MIN(84, qlen - (query_ptr - query)));
					query_ptr += MIN(84, qlen - (query_ptr - query));
				}
				memcpy(target_ptr, &buf[84], MIN(84, tlen - (target_ptr - target)));
				target_ptr += MIN(84, tlen - (target_ptr - target));
				if (line->ty0.preamble == PACKET_END)
					break;
			}
		}

		int lscore, gscore, tle, qle, gtle, max_off;
		lscore = ksw_extend2(qlen, query, tlen, target, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, opt->pen_clip5, opt->zdrop, init_score, &qle, &tle, &gtle, &gscore, &max_off);

		struct ResultLine *rl;
		if (results.size() == 0 || results.back().ty_r.preamble[0] >= 5)
		{
			results.push_back({0});
			rl = &results.back().ty_r;
		}
		rl = &results.back().ty_r;

		struct ResultEntry &re = rl->results[rl->preamble[0]++];
		re.qle = qle;
		re.gscore = gscore;
		re.gtle = gtle;
		re.lscore = lscore;
		re.seq_id = seq_id;
		re.tle = tle;
		free(query);
		free(target);
		// test exception
		re.exception = 0;
	}
}

void dump_mem(const char *fname, const std::vector<union SeedExLine>& buf)
{
	FILE *fp = fopen(fname, "w");
	assert(fname);
	fwrite(buf.data(), sizeof(union SeedExLine), buf.size(), fp);
	fclose(fp);
}

void mem_chain2aln_cpu(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av, int64_t rmax0, int64_t rmax1)
{
	int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp, max = 0;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	uint64_t *srt;

	if (c->n == 0) return;
		// FPGA : Write read data into write_buffer
		if (bwa_verbose >= 10) {
				int j;
				printf("*** FPGA : Seeing Read Query:   "); for (j = 0; j < l_query; ++j) putchar("ACGTN"[(int)query[j]]); putchar('\n');
		}


		rmax[0] = rmax0;
		rmax[1] = rmax1;

	// retrieve the reference sequence
	rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
	assert(c->rid == rid);

	// for (k = c->n - 1; k >= 0; --k) {
	for (k = 0; k < c->n; ++k) {
		mem_alnreg_t *a;
		// s = &c->seeds[(uint32_t)srt[k]];                    // Select seed with best score first within a chain
		s = &c->seeds[k];                    // Select seed with best score first within a chain

		a = kv_pushp(mem_alnreg_t, *av);
		memset(a, 0, sizeof(mem_alnreg_t));
		a->w = aw[0] = aw[1] = opt->w;
		a->score = a->truesc = -1;
		a->rid = c->rid;

		if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);
				if(bwa_verbose >= 10) {
					printf("FPGA Qbeg : %d\n",s->qbeg);
					printf("FPGA Qend : %d\n",s->qbeg + s->len);
					printf("FPGA Seed beg : %ld\n",s->rbeg);
				}

		if (s->qbeg) { // left extension
			uint8_t *rs, *qs;
			int qle, tle, gtle, gscore;
			qs = malloc(s->qbeg);
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = malloc(tmp);
			for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[0] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
					printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
				}
				a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
				if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
				if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
			}
			// check whether we prefer to reach the end of the query
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
				a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
				a->truesc = a->score;
			} else { // to-end extension
				a->qb = 0, a->rb = s->rbeg - gtle;
				a->truesc = gscore;
			}
			free(qs); free(rs);
		} else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

		if (s->qbeg + s->len != l_query) { // right extension
			int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
			qe = s->qbeg + s->len;
			re = s->rbeg + s->len - rmax[0];
			assert(re >= 0);
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[1] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
					printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
				}
				a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
				if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
				if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
			}
			// similar to the above
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
				a->qe = qe + qle, a->re = rmax[0] + re + tle;
				a->truesc += a->score - sc0;
			} else { // to-end extension
				a->qe = l_query, a->re = rmax[0] + re + gtle;
				a->truesc += gscore - sc0;
			}
		} else a->qe = l_query, a->re = s->rbeg + s->len;
		if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);

		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
		a->w = aw[0] > aw[1]? aw[0] : aw[1];
		a->seedlen0 = s->len;

		a->frac_rep = c->frac_rep;
	}
	free(rseq);
}


//void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av)
void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av, int64_t rmax0, int64_t rmax1)
{
	int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp, max = 0;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	uint64_t *srt;


	if (c->n == 0) return;
		// FPGA : Write read data into write_buffer
		if (bwa_verbose >= 10) {
				int j;
				printf("*** FPGA : Seeing Read Query:   "); for (j = 0; j < l_query; ++j) putchar("ACGTN"[(int)query[j]]); putchar('\n');
		}
		
			 
	// get the max possible span
	/*rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len;
	}
	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;


	if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
		if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else rmax[0] = l_pac;
	}
		if (bwa_verbose >= 10) {
				int j;
				printf("*** FPGA : rmax[0] :%llu, rmax[1]: %llu   \n",rmax[0],rmax[1]);
		}*/
		
		rmax[0] = rmax0;
		rmax[1] = rmax1;

	// retrieve the reference sequence
	rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
	assert(c->rid == rid);

	srt = malloc(c->n * 8);
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].score<<32 | i;       // Fill srt with first 32 bits as the index in chain , and the higher 32 bits as score for sorting
	ks_introsort_64(c->n, srt);

	for (k = c->n - 1; k >= 0; --k) {
		mem_alnreg_t *a;
		s = &c->seeds[(uint32_t)srt[k]];                    // Select seed with best score first within a chain

		for (i = 0; i < av->n; ++i) { // test whether extension has been made before
			mem_alnreg_t *p = &av->a[i];
			int64_t rd;
			int qd, w, max_gap;
			if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
			if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
			// qd: distance ahead of the seed on query; rd: on reference
			qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
			max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
			w = max_gap < p->w? max_gap : p->w; // bounded by the band width
			if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
			// similar to the previous four lines, but this time we look at the region behind
			qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
			max_gap = cal_max_gap(opt, qd < rd? qd : rd);
			w = max_gap < p->w? max_gap : p->w;
			if (qd - rd < w && rd - qd < w) break;
		}
		if(bwa_verbose >= 18){
				printf("(FPGA) i = %d,k = %d, av_size = %zu\n",i,k,av->n);
		}
		if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
			if (bwa_verbose >= 4)
				printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
					   k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
			for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
				const mem_seed_t *t;
				if (srt[i] == 0) continue;
				t = &c->seeds[(uint32_t)srt[i]];
				if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
				if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
				if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
			}
			if (i == c->n) { // no overlapping seeds; then skip extension
				srt[k] = 0; // mark that seed extension has not been performed
				continue;
			}
			if (bwa_verbose >= 4)
				printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
		}

		a = kv_pushp(mem_alnreg_t, *av);
		memset(a, 0, sizeof(mem_alnreg_t));
		a->w = aw[0] = aw[1] = opt->w;
		a->score = a->truesc = -1;
		a->rid = c->rid;

		if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);
				if(bwa_verbose >= 10) {
					printf("FPGA Qbeg : %d\n",s->qbeg);
					printf("FPGA Qend : %d\n",s->qbeg + s->len);
					printf("FPGA Seed beg : %ld\n",s->rbeg);
				}

		if (s->qbeg) { // left extension
				


			uint8_t *rs, *qs;
			int qle, tle, gtle, gscore;
			qs = malloc(s->qbeg);
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = malloc(tmp);
			for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[0] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
					printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
				}
				a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
				if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
				if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
			}
			// check whether we prefer to reach the end of the query
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
				a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
				a->truesc = a->score;
			} else { // to-end extension
				a->qb = 0, a->rb = s->rbeg - gtle;
				a->truesc = gscore;
			}
			free(qs); free(rs);
		} else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

		if (s->qbeg + s->len != l_query) { // right extension
			int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
			qe = s->qbeg + s->len;
			re = s->rbeg + s->len - rmax[0];
			assert(re >= 0);
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[1] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
					printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
				}
				a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
				if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
				if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
			}
			// similar to the above
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
				a->qe = qe + qle, a->re = rmax[0] + re + tle;
				a->truesc += a->score - sc0;
			} else { // to-end extension
				a->qe = l_query, a->re = rmax[0] + re + gtle;
				a->truesc += gscore - sc0;
			}
		} else a->qe = l_query, a->re = s->rbeg + s->len;
		if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);

		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
		a->w = aw[0] > aw[1]? aw[0] : aw[1];
		a->seedlen0 = s->len;

		a->frac_rep = c->frac_rep;
	}
	free(srt); free(rseq);
}

/*****************************
 * Basic hit->SAM conversion *
 *****************************/

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r)
{
	int w;
	if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
	w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
	if (w < abs(l1 - l2)) w = abs(l1 - l2);
	return w;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k) {
		int op = cigar[k]&0xf;
		if (op == 0 || op == 2)
			l += cigar[k]>>4;
	}
	return l;
}

static inline void add_cigar(const mem_opt_t *opt, mem_aln_t *p, kstring_t *str, int which)
{
	int i;
	if (p->n_cigar) { // aligned
		for (i = 0; i < p->n_cigar; ++i) {
			int c = p->cigar[i]&0xf;
			if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
				c = which? 4 : 3; // use hard clipping for supplementary alignments
			kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
		}
	} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
}

//Only support for single ended reads
/*void mem_aln2bin_sam(const mem_opt_t *opt, const bntseq_t *bns, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_){
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	s->abs_pos[which] = p->aln_abs_pos;
	s->correction[which] = p->correction;
	s->is_rev[which] = p->is_rev;
}*/

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, int64_t *abs_pos)
{
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand




	// print up to CIGAR
	l_name = strlen(s->name);
	ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
	kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		add_cigar(opt, p, str, which);
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);

	// put abs position in the sequence record
	s->abs_pos = p->aln_abs_pos;
	s->correction = p->correction;
	s->is_rev = p->is_rev;


	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}

	// print optional tags
	if (p->n_cigar) {
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (m && m->n_cigar) { kputsn("\tMC:Z:", 6, str); add_cigar(opt, m, str, which); }
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str);
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
		if (p->alt_sc > 0)
			ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
	}
	if (p->XA) {
		kputsn((opt->flag&MEM_F_XB)? "\tXB:Z:" : "\tXA:Z:", 6, str);
		kputs(p->XA, str);
	}
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
		int tmp;
		kputsn("\tXR:Z:", 6, str);
		tmp = str->l;
		kputs(bns->anns[p->rid].anno, str);
		for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
			if (str->s[i] == '\t') str->s[i] = ' ';
	}
	kputc('\n', str);
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a)
{
	int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
	double identity;
	sub = a->csub > sub? a->csub : sub;
	if (sub >= a->score) return 0;
	l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
	identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
	if (a->score == 0) {
		mapq = 0;
	} else if (opt->mapQ_coef_len > 0) {
		double tmp;
		tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
		tmp *= identity * identity;
		mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
	} else {
		mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
		mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	}
	if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
	if (mapq > 60) mapq = 60;
	if (mapq < 0) mapq = 0;
	mapq = (int)(mapq * (1. - a->frac_rep) + .499);
	return mapq;
}

void mem_reorder_primary5(int T, mem_alnreg_v *a)
{
	int k, n_pri = 0, left_st = INT_MAX, left_k = -1;
	mem_alnreg_t t;
	for (k = 0; k < a->n; ++k)
		if (a->a[k].secondary < 0 && !a->a[k].is_alt && a->a[k].score >= T) ++n_pri;
	if (n_pri <= 1) return; // only one alignment
	for (k = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		if (p->secondary >= 0 || p->is_alt || p->score < T) continue;
		if (p->qb < left_st) left_st = p->qb, left_k = k;
	}
	assert(a->a[0].secondary < 0);
	if (left_k == 0) return; // no need to reorder
	t = a->a[0], a->a[0] = a->a[left_k], a->a[left_k] = t;
	for (k = 1; k < a->n; ++k) { // update secondary and secondary_all
		mem_alnreg_t *p = &a->a[k];
		if (p->secondary == 0) p->secondary = left_k;
		else if (p->secondary == left_k) p->secondary = 0;
		if (p->secondary_all == 0) p->secondary_all = left_k;
		else if (p->secondary_all == left_k) p->secondary_all = 0;
	}
}


/*void create_bin_sam(bseq1_t * s){
	int l_seq = s->l_seq;
	int i = 0;
	s->bin_sam = (uint8_t *) malloc(l_seq * sizeof(uint8_t));
	memset(s->bin_sam,0,(l_seq * sizeof(uint8_t)));
	for(i = 0;i<l_seq;i++){
		if(s->seq[i] < 4){
			s->bin_sam[i] = ((int)s->qual[i] - 33) << 2;
			s->bin_sam[i] |= s->seq[i];
		}
	}
}*/


// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
	
	extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, mem_alnreg_v *a, int l_query, const char *query);
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k, l;
	char **XA = 0;
	

	int64_t abs_pos = 0;
	if (!(opt->flag & MEM_F_ALL))
		XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
	kv_init(aa);
	str.l = str.m = 0; str.s = 0;
	for (k = l = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;
		if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		assert(q->rid >= 0); // this should not happen with the new code
		q->XA = XA? XA[k] : 0;
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (l && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && l && !p->is_alt && q->mapq > aa.a[0].mapq)
			q->mapq = aa.a[0].mapq; // lower mapq for supplementary mappings, unless -5 or -q is applied
		++l;
	}


	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
			mem_aln_t t;
			t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
			t.flag |= extra_flag;
			//mem_aln2bin_sam(opt, bns, s, 1, &t, 0, m);
			mem_aln2sam(opt, bns, &str, s, 1, &t, 0, m, &abs_pos);
			//s->abs_pos = abs_pos;
	} else {
			for (k = 0; k < aa.n; ++k){
				//mem_aln2bin_sam(opt, bns, s, aa.n, aa.a, k, m);
				mem_aln2sam(opt, bns, &str, s, aa.n, aa.a, k, m, &abs_pos);
			}
			for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
			free(aa.a);
			//s->abs_pos = abs_pos;
	}
	s->sam = str.s;
	if (XA) {
			for (k = 0; k < a->n; ++k) free(XA[k]);
			free(XA);
	}
}



void seeding(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf, mem_chain_v ** chain){
	int i;
	mem_chain_v * chn = (mem_chain_v *) malloc(sizeof(mem_chain_v));


	for (i = 0; i < l_seq; ++i) // convert to 2-bit encoding if we have not done so
		seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]];

	*chn = mem_chain(opt, bwt, bns, l_seq, (uint8_t*)seq, buf);
	chn->n = mem_chain_flt(opt, chn->n, chn->a);
	mem_flt_chained_seeds(opt, bns, pac, l_seq, (uint8_t*)seq, chn->n, chn->a);
	if (bwa_verbose >= 4) mem_print_chain(bns, chn);

	*chain = chn;

	return;
}


void seed_extension(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, mem_chain_v * chn, mem_alnreg_v_v * alnregs, fpga_data_out_t *data_out, std::vector<union SeedExLine>& write_buffer1, std::vector<union SeedExLine*>& write_buffer_entry_idx1, std::vector<union SeedExLine>& write_buffer2, std::vector<union SeedExLine*>& write_buffer_entry_idx2, std::vector<struct extension_meta_t>& extension_meta, int run_fpga){
	
	// mem_alnreg_v * regs = (mem_alnreg_v *) malloc(sizeof(mem_alnreg_v));
	// memset(regs,0,sizeof(mem_alnreg_v));
	// kv_init(*alnregs);
	mem_alnreg_v * regs = NULL;

	int i = 0;
	for (i = 0; i < chn->n; ++i) {
		mem_chain_t *p = &(chn->a[i]);
		if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
		int64_t rmax0 = 0;
		int64_t rmax1 = 0;

		fetch_rmaxs(opt, bns,pac, l_seq, (uint8_t*) seq, p, regs, &rmax0, &rmax1);
		if(run_fpga == 0){
			if (!regs) {
				if (alnregs) {
					for (int j = 0; j < alnregs->n; j++) free(alnregs->a[j].a);
					// kv_resize(mem_alnreg_v, *alnregs, 0);
					alnregs->n = 0;
				}
				regs = kv_pushp(mem_alnreg_v, *alnregs);
				memset(regs, 0, sizeof(mem_alnreg_v));
			}
			mem_chain2aln(opt, bns, pac, l_seq, (uint8_t*)seq, p, regs,rmax0,rmax1);
		}
		else {
			regs = kv_pushp(mem_alnreg_v, *alnregs);
			memset(regs, 0, sizeof(mem_alnreg_v));
			if ((rmax1 - rmax0) > 250 || run_fpga == 2 /* for debug */) {
  			  	// Max allowed ref length
				mem_chain2aln_cpu(opt, bns, pac, l_seq, (uint8_t*)seq, p, regs,rmax0,rmax1);
			}
  		  	else{
  	 			extension_meta.back().chain_id = i;
  	 			mem_chain2aln_to_fpga(opt, bns, pac, l_seq, (uint8_t*)seq, p, regs,rmax0,rmax1, data_out, write_buffer1, write_buffer_entry_idx1, write_buffer2, write_buffer_entry_idx2, extension_meta);
  		  	}
		}
		//free(chn->a[i].seeds);
	}
	//free(chn->a);
	// regs->n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t*)seq, regs->n, regs->a);

	// if (bwa_verbose >= 4) {
	// 	err_printf("* %ld chains remain after removing duplicated chains\n", regs->n);
	// 	for (i = 0; i < regs->n; ++i) {
	// 		mem_alnreg_t *p = &(regs->a[i]);
	// 		printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
	// 	}
	// }
	// for (i = 0; i < regs->n; ++i) {
	// 	mem_alnreg_t *p = &(regs->a[i]);
	// 	if (p->rid >= 0 && bns->anns[p->rid].is_alt)
	// 		p->is_alt = 1;
	// }

	// if(data_out->fpga_entry_present == 1){
	//     if (bwa_verbose >= 10) {
	//         printf("Writing entry num for read\n");
	//     }
	// }
	// *alnreg = regs; 
	return;
}

mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf, fpga_data_out_t *data_out, uint8_t** write_buffer, size_t *write_buffer_size, uint32_t *write_buffer_index,uint32_t read_id)
{

		// Enter this method once for every read
	int i;
	mem_chain_v chn;
	mem_alnreg_v regs;

	for (i = 0; i < l_seq; ++i) // convert to 2-bit encoding if we have not done so
		seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]];

	chn = mem_chain(opt, bwt, bns, l_seq, (uint8_t*)seq, buf);
	chn.n = mem_chain_flt(opt, chn.n, chn.a);
	mem_flt_chained_seeds(opt, bns, pac, l_seq, (uint8_t*)seq, chn.n, chn.a);
	if (bwa_verbose >= 4) mem_print_chain(bns, &chn);

	kv_init(regs);
	for (i = 0; i < chn.n; ++i) {
		mem_chain_t *p = &chn.a[i];
		if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
				int64_t rmax0 = 0;
				int64_t rmax1 = 0;

				fetch_rmaxs(opt, bns,pac, l_seq, (uint8_t*) seq, p, &regs, &rmax0, &rmax1);
				if((rmax1 - rmax0) > 0){
					// Max allowed ref length
					mem_chain2aln(opt, bns, pac, l_seq, (uint8_t*)seq, p, &regs,rmax0,rmax1);
				}
				else{
					// mem_chain2aln_to_fpga(opt, bns, pac, l_seq, (uint8_t*)seq, p, &regs,rmax0,rmax1, data_out, write_buffer, write_buffer_size,write_buffer_index,&ar_index);
				}
		free(chn.a[i].seeds);
	}
	free(chn.a);
	regs.n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t*)seq, regs.n, regs.a);
	if (bwa_verbose >= 4) {
		err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
		for (i = 0; i < regs.n; ++i) {
			mem_alnreg_t *p = &regs.a[i];
			printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
		}
	}
	for (i = 0; i < regs.n; ++i) {
		mem_alnreg_t *p = &regs.a[i];
		if (p->rid >= 0 && bns->anns[p->rid].is_alt)
			p->is_alt = 1;
	}
		if(data_out->fpga_entry_present == 1){
		if (bwa_verbose >= 10) {
				printf("Writing entry num for read\n");
			}
		}
	return regs;
}

mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar)
{
	mem_aln_t a;
	int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1<<30), l_MD;
	int64_t pos, rb, re;
	uint8_t *query;

	memset(&a, 0, sizeof(mem_aln_t));
	if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
		a.rid = -1; a.pos = -1; a.flag |= 0x4;
		a.aln_abs_pos = bns->l_pac;
		return a;
	}
	qb = ar->qb, qe = ar->qe;
	rb = ar->rb, re = ar->re;
	query = malloc(l_query);
	for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
		query[i] = query_[i] < 5? query_[i] : nst_nt4_table[(int)query_[i]];
	a.mapq = ar->secondary < 0? mem_approx_mapq_se(opt, ar) : 0;
	if (ar->secondary >= 0) a.flag |= 0x100; // secondary alignment
	tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
	w2  = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
	w2 = w2 > tmp? w2 : tmp;
	if (bwa_verbose >= 4) printf("* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
	if (w2 > opt->w) w2 = w2 < ar->w? w2 : ar->w;
	i = 0; a.cigar = 0;
	do {
		free(a.cigar);
		w2 = w2 < opt->w<<2? w2 : opt->w<<2;

		a.cigar = bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*)&query[qb], rb, re, &score, &a.n_cigar, &NM);
		if (bwa_verbose >= 4) printf("* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
		if (score == last_sc || w2 == opt->w<<2) break; // it is possible that global alignment and local alignment give different scores
		last_sc = score;
		w2 <<= 1;
	} while (++i < 3 && score < ar->truesc - opt->a);
	l_MD = strlen((char*)(a.cigar + a.n_cigar)) + 1;
	a.NM = NM;
	pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);

	a.is_rev = is_rev;
	if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
		if ((a.cigar[0]&0xf) == 2) {
			pos += a.cigar[0]>>4;
			--a.n_cigar;
			memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
		} else if ((a.cigar[a.n_cigar-1]&0xf) == 2) {
			--a.n_cigar;
			memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
		}
	}
	int64_t inter_pos = pos;
	if (qb != 0 || qe != l_query) { // add clipping to CIGAR
		int clip5, clip3;
		clip5 = is_rev? l_query - qe : qb;
		clip3 = is_rev? qb : l_query - qe;
		inter_pos -= clip5;
		a.cigar = realloc(a.cigar, 4 * (a.n_cigar + 2) + l_MD);
		if (clip5) {
			memmove(a.cigar+1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
			a.cigar[0] = clip5<<4 | 3;
			++a.n_cigar;
		}
		if (clip3) {
			memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
			a.cigar[a.n_cigar++] = clip3<<4 | 3;
		}
	}

	a.aln_abs_pos = inter_pos;             // Get abs position of aln for sorting
	a.correction = pos - inter_pos;
	a.rid = bns_pos2rid(bns, pos);
	//a.aln_abs_pos = pos;             // Get abs position of aln for sorting
	//assert(a.rid == ar->rid);
	a.pos = pos - bns->anns[a.rid].offset;
	a.score = ar->score; a.sub = ar->sub > ar->csub? ar->sub : ar->csub;
	a.is_alt = ar->is_alt; a.alt_sc = ar->alt_sc;
	free(query);
	return a;
}



void get_scores_left(const mem_opt_t *opt, ResultEntry *re,const bntseq_t *bns, const mem_chain_t *c, mem_alnreg_v *in_a, uint32_t reg_id, bool *need_rerun){

	if (re->exception & 0x5) {
		*need_rerun = true;
		return;
	}

	// mem_alnreg_t *a = kv_pushp(mem_alnreg_t, *in_a);
	// memset(a, 0, sizeof(mem_alnreg_t));
	mem_alnreg_t *a = &in_a->a[reg_id];

	const mem_seed_t *s = &c->seeds[(uint32_t)reg_id];

	//in_fpga_result_entry->read_id = read_buffer[0] & mask; 
	if(bwa_verbose >= 15){
		printf("Read ID : %x\n",re->seq_id);
	}

	int lscore = (~re->lscore == 0)? -1 : re->lscore;
	int gscore = (~re->gscore == 0)? -1 : re->gscore;

	// store returned value
	a->score = lscore;

	// check whether we prefer to reach the end of the query
	if (gscore <= 0 || gscore <= a->score - opt->pen_clip5)
	{ // local extension
		a->qb = s->qbeg - re->qle, a->rb = s->rbeg - re->tle;
		a->truesc = a->score;
	}
	else
	{ // to-end extension
		a->qb = 0, a->rb = s->rbeg - re->gtle;
		a->truesc = gscore;
	}

}

void get_scores_right(const mem_opt_t *opt, ResultEntry *re,const bntseq_t *bns, int l_query, const mem_chain_t *c, mem_alnreg_v *in_a, uint32_t reg_id, bool *need_rerun){

	if (re->exception & 0x5) {
		*need_rerun = true;
		return;
	}

	mem_alnreg_t *a = &in_a->a[reg_id];

	const mem_seed_t *s = &c->seeds[(uint32_t)reg_id];
	int sc0 = a->score;

	//in_fpga_result_entry->read_id = read_buffer[0] & mask; 
	if(bwa_verbose >= 15){
		printf("Read ID : %x\n",re->seq_id);
	}

	int lscore = (~re->lscore == 0)? -1 : re->lscore;
	int gscore = (~re->gscore == 0)? -1 : re->gscore;

	// store returned value
	a->score = lscore;

	// similar to the above
	if (gscore <= 0 || gscore <= a->score - opt->pen_clip3)
	{ // local extension
		a->qe = s->qbeg + s->len + re->qle, a->re = s->rbeg + s->len + re->tle;
		a->truesc += a->score - sc0;
	}
	else
	{ // to-end extension
		a->qe = l_query, a->re = s->rbeg + s->len + re->gtle;
		a->truesc += gscore - sc0;
	}

	// compute seedcov
	int i;
	for (i = 0, a->seedcov = 0; i < c->n; ++i) {
		const mem_seed_t *t = &c->seeds[i];
		if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
			a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
	}
}

void rerun_left_extension(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, mem_chain_t *c, mem_alnreg_v *av, int reg_id){
	mem_alnreg_t *a = &av->a[reg_id];
	const mem_seed_t *s = &c->seeds[(uint32_t)reg_id];
	if (s->qbeg) { // left extension
		// printf("@@@ Rerunning rerun left\n");
		int i, rid, aw[2], max_off[2];
		int64_t rmax[2], tmp;
		fetch_rmaxs(opt, bns,pac, l_query, (uint8_t*) query, c, NULL, &rmax[0], &rmax[1]);
		uint8_t *rseq = 0;
		rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
		uint8_t *rs, *qs;
		int qle, tle, gtle, gscore;
		qs = malloc(s->qbeg);
		for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
		tmp = s->rbeg - rmax[0];
		rs = malloc(tmp);
		for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
		for (i = 0; i < MAX_BAND_TRY; ++i) {
			int prev = a->score;
			aw[0] = opt->w << i;
			if (bwa_verbose >= 4) {
				int j;
				printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
				printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
			}
			a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
			if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
			if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
		}
		// check whether we prefer to reach the end of the query
		if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
			a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
			a->truesc = a->score;
		} else { // to-end extension
			a->qb = 0, a->rb = s->rbeg - gtle;
			a->truesc = gscore;
		}
		free(qs); free(rs);
	} else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;
}

void rerun_right_extension(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, mem_chain_t *c, mem_alnreg_v *av, int reg_id){
	mem_alnreg_t *a = &av->a[reg_id];
	const mem_seed_t *s = &c->seeds[(uint32_t)reg_id];
	if (s->qbeg + s->len != l_query) { // right extension
		// printf("@@@ Rerunning rerun right\n");
		int i, rid, aw[2], max_off[2];
		int64_t rmax[2];
		fetch_rmaxs(opt, bns,pac, l_query, (uint8_t*) query, c, NULL, &rmax[0], &rmax[1]);
		uint8_t *rseq = 0;
		rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
		int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
		qe = s->qbeg + s->len;
		re = s->rbeg + s->len - rmax[0];
		assert(re >= 0);
		for (i = 0; i < MAX_BAND_TRY; ++i) {
			int prev = a->score;
			aw[1] = opt->w << i;
			if (bwa_verbose >= 4) {
				int j;
				printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
				printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
			}
			a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
			if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
			if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
		}
		// similar to the above
		if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
			a->qe = qe + qle, a->re = rmax[0] + re + tle;
			a->truesc += a->score - sc0;
		} else { // to-end extension
			a->qe = l_query, a->re = rmax[0] + re + gtle;
			a->truesc += gscore - sc0;
		}
	} else a->qe = l_query, a->re = s->rbeg + s->len;
}

void get_all_scores(const worker_t *w, uint8_t *read_buffer, int total_lines, queue_t *qe,fpga_data_out_v * f1v, std::vector<struct extension_meta_t>& extension_meta, mem_alnreg_v_v *alnregs){
	int i = 0;
	for(i=0;i<total_lines;i++){
		if(bwa_verbose >= 15){
			int k1=0,i1=0;
			for(k1 = 63;k1>=0;k1--) {
				for(i1 = 8-4;i1>=0;i1 -= 4){
					printf("%x",read_buffer[i*64 + k1]>>i1 & 0xF);
				}

			}
			printf("\n");
		}

		static_assert(sizeof(ResultLine) == sizeof(SeedExLine), "Line size mismatch");
		struct ResultLine* results = ((struct ResultLine*) read_buffer) + i;

		// memcpy(&read_id,read_buffer + i*64 + 1,4);
		// read_id = read_id - qe->starting_read_id;

		for (int k = 0; k < results->preamble[0]; ++k) {
			struct ResultEntry * re = &results->results[k];
			int seq_id = re->seq_id;
			uint32_t read_idx = extension_meta.at(re->seq_id).read_idx;
			// uint32_t read_id = qe->seqs[read_idx]->read_id;
			uint32_t chain_id = extension_meta.at(re->seq_id).chain_id;
			uint32_t seed_id = extension_meta.at(re->seq_id).seed_id;
			//TODO: Fix this if condition
			assert(read_idx < BATCH_SIZE);
			assert(f1v->a[read_idx].fpga_entry_present == 1);
			// if(read_idx < BATCH_SIZE){
				// if(f1v->a[read_idx].fpga_entry_present == 1){
					bool need_rerun = false;
					if (f1v->read_right) {
						get_scores_right(w->opt, re, global_bns, qe->seqs[read_idx]->l_seq,&qe->chains[read_idx]->a[chain_id], &(alnregs[read_idx].a[chain_id]), seed_id, &need_rerun);
						if (need_rerun) {
							rerun_right_extension(w->opt, w->bwt, w->bns, w->pac, qe->seqs[read_idx]->l_seq, (const uint8_t *)qe->seqs[read_idx]->seq, &qe->chains[read_idx]->a[chain_id], &(alnregs[read_idx].a[chain_id]), seed_id);
						}
					} else {
						get_scores_left(w->opt, re, global_bns,&qe->chains[read_idx]->a[chain_id], &(alnregs[read_idx].a[chain_id]), seed_id, &need_rerun);
						if (need_rerun) {
							rerun_left_extension(w->opt, w->bwt, w->bns, w->pac, qe->seqs[read_idx]->l_seq, (const uint8_t *)qe->seqs[read_idx]->seq, &qe->chains[read_idx]->a[chain_id], &(alnregs[read_idx].a[chain_id]), seed_id);
						}

						// transfer score for sc0 in loadbuf2
						union SeedExLine * right_ext_entry;
						if (right_ext_entry = f1v->load_buffer_entry_idx2->at(re->seq_id)) {
							right_ext_entry->ty1.params.init_score = re->lscore;
						}
					}
				// }
			// }
		}

	}

}


void read_scores_from_fpga(const worker_t *w, pci_bar_handle_t pci_bar_handle,queue_t* qe, fpga_data_out_v * f1v, int channel, uint64_t addr, std::vector<struct extension_meta_t>& extension_meta, mem_alnreg_v_v *alnregs){
	int rc = 0;
	 
	if(f1v->n != 0){   
		if(bwa_verbose >= 10) {
			printf("Num entries in read_from_fpga : %zd\n",f1v->n);
		}

		size_t total_lines = ((f1v->read_right? f1v->load_buffer_entry_idx2->size() : f1v->load_buffer_entry_idx1->size()) - 2 + (sizeof(ResultLine::results) / sizeof(ResultEntry))) / (sizeof(ResultLine::results) / sizeof(ResultEntry)) ;
		size_t read_buffer_size = total_lines * 64;

#ifdef ENABLE_FPGA
		uint8_t * read_buffer = read_from_fpga(fpga_pci_local->read_fd,read_buffer_size,channel * MEM_16G + addr);

		get_all_scores(w,read_buffer,total_lines,qe,f1v,extension_meta, alnregs);

		if(read_buffer) {
			free(read_buffer);
		}
#endif
	}
	
}

void delete_queue_entry(queue_t *qe){
	if(qe == NULL){
		return;
	}


	// Pass on qe to the next stage

	if(qe->chains) {
		free(qe->chains);
	}
	if(qe->regs) {
		int i = 0;
		for(i=0;i<qe->num;i++){
			if(qe->regs[i]->a)
				free(qe->regs[i]->a);
		}
		free(qe->regs);
	}
	if(qe->seqs){
		free(qe->seqs);
	}
	free(qe);
	qe = NULL;
	return;
}












/*typedef struct {
	worker_t *w_master;     // Location of master with all sequences
	queue *queue1;          // Location of common queue for fpga
	int tid;                // Thread id

	// Sequences processed by any thread will be all seqs starting from tid*BATCH_SIZE;
	// next batch to be processed will be opt->n_threads*BATCH_SIZE

} worker_slave_t;*/








void worker1_ST(void *data){

	worker_slave_t *slave_data = (worker_slave_t*)data;
	int tid = slave_data->tid;

	worker_t *w = slave_data->w_master;
	queue *q = w->queue1;
	int n_threads = w->opt->n_threads;

	int64_t i = 0;
	int j = 0;

	queue_t *qe;


	/*if(w->n_processed == 0){
	}*/

	// one thread will process k=(tot+nth-1)/nth seeds indexed as [tid * k, (tid + 1)k)
	int K = (w->n_processed + n_threads - 1) / n_threads;
	for(i = tid*K; i < (tid + 1) * K; i += j){
		if (i>=w->n_processed) break;

		qe = (queue_t*)malloc(sizeof(queue_t));
		qe->regs = (mem_alnreg_v **)malloc(BATCH_SIZE * sizeof(mem_alnreg_v *));
		qe->chains = (mem_chain_v **)malloc(BATCH_SIZE * sizeof(mem_chain_v *));
		qe->seqs = (bseq1_t **)malloc(BATCH_SIZE * sizeof(bseq1_t *));
		qe->num = 0;
		qe->last_entry = 0;
		qe->starting_read_id = i;
		for(j = 0;j<BATCH_SIZE;j++){
			int n_lines = 0;

			if(i+j < w->n_processed){
				w->seqs[i+j].read_id = i+j;
				qe->seqs[j] = &w->seqs[i+j];
				qe->num++;

				if (!(w->opt->flag&MEM_F_PE)) {
						if (bwa_verbose >= 4) printf("=====> Processing read '%s'| (i+j) = %ld  <=====\n", w->seqs[i+j].name,(i+j));
						
						seeding(w->opt, w->bwt, w->bns, w->pac, w->seqs[i+j].l_seq, w->seqs[i+j].seq, w->aux[tid], &qe->chains[j]);
						n_lines += qe->chains[j]->n * 3;
//                        qe->regs[j] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i+j].l_seq, w->seqs[i+j].seq, w->aux[tid],&f1, &qe->load_buffer, &(qe->load_buffer_size),&write_buffer_index,(uint32_t)(i+j));
						//qe->fpga_results->a[j].fpga_entry_present = f1.fpga_entry_present;
						/*if(f1.fpga_entry_present){
							qe->fpga_results->n++;
						}*/
						//smem_aux_destroy(aux1);

				} else {
						if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1 <=====\n", w->seqs[(i+j)<<1|0].name);
						seeding(w->opt, w->bwt, w->bns, w->pac, w->seqs[(i+j)<<1|0].l_seq, w->seqs[(i+j)<<1|0].seq, w->aux[tid], &qe->chains[j<<1|0]);
						n_lines += qe->chains[j<<1|0]->n * 3;
						//qe->regs[(j)<<1|0] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[(i+j)<<1|0].l_seq, w->seqs[(i+j)<<1|0].seq, w->aux[tid],&f1,&qe->load_buffer, &(qe->load_buffer_size),&write_buffer_index,(uint32_t)((i+j)<<1|0));
						//qe->fpga_results->a[j<<1|0].fpga_entry_present = f1.fpga_entry_present;
						//if(f1.fpga_entry_present){
						//    qe->fpga_results->n++;
						//}
						//smem_aux_destroy(aux1);
						
						if (bwa_verbose >= 4) printf("=====> Processing read '%s'/2 <=====\n", w->seqs[(i+j)<<1|1].name);

						seeding(w->opt, w->bwt, w->bns, w->pac, w->seqs[(i+j)<<1|1].l_seq, w->seqs[(i+j)<<1|1].seq, w->aux[tid], &qe->chains[j<<1|1]);
						n_lines += qe->chains[j<<1|1]->n * 3;
						//qe->regs[(j)<<1|1] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[(i+j)<<1|1].l_seq, w->seqs[(i+j)<<1|1].seq, w->aux[tid],&f1,&qe->load_buffer, &(qe->load_buffer_size),&write_buffer_index,(uint32_t)((i+j)<<1|1));
						//qe->fpga_results->a[i<<1|1].fpga_entry_present = f1.fpga_entry_present;
						//if(f1.fpga_entry_present){
						//    qe->fpga_results->n++;
						//}
						//smem_aux_destroy(aux1);
				}

				if (n_lines >= BATCH_LINE_LIMIT) {
					// err_printf("@@@ Limit batchsize to avoid buffer overflow (at:%d %d)\n", n_lines, j);
					assert(j > 0 && "Batch line size is too small (cannot pack even 1 read).");
					// revoke last entry
					qe->num--; j--;
					break;
				}

			}
			else{
				break;
			}
			 
		
		}
		/*if((i+BATCH_SIZE) >= w->n_processed){
			qe->last_entry = 1;
		}*/

		
		// Grab queue mutex and add queue_element in the queue
		pthread_mutex_lock (q->mut);
		while (q->full) {
			if(bwa_verbose >= 18){
				printf ("producer: queue FULL.\n");
			}
			pthread_cond_wait (q->notFull, q->mut);
		}
		queueAdd (q, qe);
		pthread_mutex_unlock (q->mut);
		pthread_cond_signal (q->notEmpty);

	}

	pthread_exit(0);
	//return;
}

void worker1_MT(void *data){
	worker_t *w = (worker_t*)data;
	queue *q = w->queue1;
	int i = 0;

	pthread_t *w1_slaves = (pthread_t*)malloc(w->opt->n_threads * sizeof(pthread_t));
	worker_slave_t **slaves = (worker_slave_t**)malloc(w->opt->n_threads * sizeof(worker_slave_t*));


	for(i = 0;i<w->opt->n_threads;i++){
		slaves[i] = (worker_slave_t*)malloc(sizeof(worker_slave_t));
		slaves[i]->w_master = w;
		slaves[i]->tid = i;
		pthread_create (&w1_slaves[i], NULL, worker1_ST, slaves[i]);
	}

	for(i = 0;i<w->opt->n_threads;i++){
		pthread_join (w1_slaves[i], NULL);
		free(slaves[i]);
	}
	
	free(slaves);
	free(w1_slaves);

	queue_t *qe;
	qe = (queue_t*)malloc(sizeof(queue_t));
	qe->num = 0;
	qe->last_entry = 1;
	qe->regs = NULL;
	qe->chains = NULL;
	qe->seqs = NULL;
	pthread_mutex_lock (q->mut);
	while (q->full) {
		if(bwa_verbose >= 18){
			printf ("producer: queue FULL.\n");
		}
		pthread_cond_wait (q->notFull, q->mut);
	}
	for (int j = 0; j < NUM_FPGA_THREADS; ++j) queueAdd (q, qe);
	pthread_mutex_unlock (q->mut);
	pthread_cond_signal (q->notEmpty);

	
	//return;
	pthread_exit(0);
}

void free_chains(mem_chain_v * chn){
	int i = 0;
	for(i=0;i<chn->n;i++){
		free(chn->a[i].seeds);
	}
	free(chn->a);
	free(chn);
}


static void fpga_worker(void *data){
	queue_coll *qc = (queue_coll *)data;
	worker_t * w = qc->w;
	queue *q1 = qc->q1;
	queue *q2 = qc->q2;
	const int tid = qc->tid;

	// Grab mutex and get head of queue

	queue_t *qe;
	int last_entry = 0;
	int rc = 0;
	std::vector<union SeedExLine> load_buffer1;
	std::vector<union SeedExLine*> load_buffer_entry_idx1;
	std::vector<union SeedExLine> load_buffer2;
	std::vector<union SeedExLine*> load_buffer_entry_idx2;
	std::vector<struct extension_meta_t> extension_meta;
	int time_out = 0;
	struct timespec start,end;
	uint64_t timediff;


	fpga_data_out_v f1v;
	fpga_data_out_t f1;

	load_buffer1.reserve(BATCH_LINE_LIMIT);
	load_buffer2.reserve(BATCH_LINE_LIMIT);
	extension_meta.reserve(BATCH_LINE_LIMIT);
	f1v.load_buffer1 = &load_buffer1;
	f1v.load_buffer_entry_idx1 = &load_buffer_entry_idx1;
	f1v.load_buffer2 = &load_buffer2;
	f1v.load_buffer_entry_idx2 = &load_buffer_entry_idx2;

	while(1){
		pthread_mutex_lock (q1->mut);
		while (q1->empty) {
			if(bwa_verbose >= 18)
				printf ("consumer: queue EMPTY.\n");
			pthread_cond_wait (q1->notEmpty, q1->mut);
		}

		queueDel (q1, &qe);
		pthread_mutex_unlock (q1->mut);
		pthread_cond_signal (q1->notFull);
		//last_entry = 0;
		last_entry = qe->last_entry;
		fpga_mem_write_offset = 0;

		if(last_entry == 0){

			time_out = 0;

			int i = 0;

			f1v.a = (fpga_data_out_t *) malloc(qe->num * sizeof(fpga_data_out_t));
			f1v.n = 0;
			extension_meta.push_back({0, 0, 0});

			mem_alnreg_v_v * alnregs = (mem_alnreg_v_v *)calloc(qe->num, sizeof(mem_alnreg_v_v)); // read->chain->reg
			// mem_alnreg_v_v * alnregs_vv = (mem_alnreg_v_v *)calloc(qe->num, sizeof(mem_alnreg_v_v)); // read->chain->reg

			for(i = 0;i<qe->num;i++){
				f1.fpga_entry_present = 0;
				extension_meta.back().read_idx = i;
				kv_init(alnregs[i]);
				seed_extension(w->opt, w->bwt, w->bns, w->pac, qe->seqs[i]->l_seq, qe->seqs[i]->seq, qe->chains[i], &alnregs[i], &f1, load_buffer1, load_buffer_entry_idx1, load_buffer2, load_buffer_entry_idx2, extension_meta, 1);
				// kv_init(alnregs_vv[i]);
				// seed_extension(w->opt, w->bwt, w->bns, w->pac, qe->seqs[i]->l_seq, qe->seqs[i]->seq, qe->chains[i], &alnregs_vv[i], &f1, load_buffer1, load_buffer_entry_idx1, load_buffer2, load_buffer_entry_idx2, extension_meta, 2);

				f1v.a[i].fpga_entry_present = f1.fpga_entry_present;
				if(f1.fpga_entry_present){
				  f1v.n++;
				}
				// Dont free chains yet
				//free(qe->chains[i]);
			}

			// push sentinel
			load_buffer1.push_back({PACKET_COMPLETE});
			load_buffer2.push_back({PACKET_COMPLETE});

			if(f1v.n != 0){
				// load_buffer = (uint8_t *)realloc(load_buffer,load_buffer_size + write_buffer_capacity);
				// memset(load_buffer + load_buffer_size,0,write_buffer_capacity);

#ifdef ENABLE_FPGA
				write_to_fpga(fpga_pci_local->write_fd,(uint8_t*)load_buffer1.data(),load_buffer1.size() * sizeof(union SeedExLine),0);

				// vdip = 0x0001;
				vdip = 2 * tid + 1;

				pthread_mutex_lock (qc->seedex_mut);

				// PCI Poke can be used for writing small amounts of data on the OCL bus
				rc = fpga_pci_poke(fpga_pci_local->pci_bar_handle,0,vdip);

				clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
				while(1) {

					rc = fpga_pci_peek(fpga_pci_local->pci_bar_handle,0,&vled);

					if(vled == 0x10)  {
						vdip = 0x0000;
						rc = fpga_pci_poke(fpga_pci_local->pci_bar_handle,0,vdip);
						break;
					}

					clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
					timediff = (end.tv_nsec - start.tv_nsec) * 1000000000 + (end.tv_sec - start.tv_sec);
					if(timediff > TIMEOUT){
						if(bwa_verbose >= 10){
							fprintf(stderr,"Going into timeout mode\n");
							fprintf(stderr,"Starting : %ld\n",qe->starting_read_id);
						}
						vdip = 0xffffffff;
						rc = fpga_pci_poke(fpga_pci_local->pci_bar_handle,0,vdip);
						do { fpga_pci_peek(fpga_pci_local->pci_bar_handle,0,&vled); } while (vled != 0x0);
						time_out = 1;
						break;
					}
				}

				pthread_mutex_unlock (qc->seedex_mut);

				if(time_out == 0){
					f1v.read_right = false;
					read_scores_from_fpga(w, bw_pci_bar_handle,qe,&f1v,0, BATCH_LINE_LIMIT*64*4 + (2*tid) * BATCH_LINE_LIMIT/4*64, extension_meta, alnregs);
				}
#else
				std::vector<union SeedExLine> read_buffer;
				f1v.read_right = false;
				pthread_mutex_lock (qc->seedex_mut);
				fpga_func_model(w->opt, load_buffer1, load_buffer_entry_idx1, read_buffer);
				pthread_mutex_unlock (qc->seedex_mut);

				// static bool dumped = false;
				// if (!dumped){
				// 	fprintf(stderr, "Dumping in.out(%d lines) out.mem(%d lines)...", load_buffer1.size(), read_buffer.size());
				// 	dump_mem("in.mem", load_buffer1);
				// 	dump_mem("out.mem", read_buffer);
				// }
				// dumped = true;

				// FIXME:: Count only non-null from entry idx
				// assert(read_buffer.size() == (load_buffer_entry_idx1.size() - 2 + (sizeof(ResultLine::results) / sizeof(ResultEntry))) / (sizeof(ResultLine::results) / sizeof(ResultEntry)) ) ;
				get_all_scores(w,(uint8_t *)read_buffer.data(),read_buffer.size(),qe,&f1v,extension_meta, alnregs);
#endif



#ifdef ENABLE_FPGA
				// right ext
				write_to_fpga(fpga_pci_local->write_fd,(uint8_t*)load_buffer2.data(),load_buffer2.size() * sizeof(union SeedExLine),0);

				// vdip = 0x0001;
				vdip = 2 * tid + 2;

				pthread_mutex_lock (qc->seedex_mut);

				// PCI Poke can be used for writing small amounts of data on the OCL bus
				rc = fpga_pci_poke(fpga_pci_local->pci_bar_handle,0,vdip);

				clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
				while(1) {

					rc = fpga_pci_peek(fpga_pci_local->pci_bar_handle,0,&vled);

					if(vled == 0x10)  {
						vdip = 0x0000;
						rc = fpga_pci_poke(fpga_pci_local->pci_bar_handle,0,vdip);
						break;
					}

					clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
					timediff = (end.tv_nsec - start.tv_nsec) * 1000000000 + (end.tv_sec - start.tv_sec);
					if(timediff > TIMEOUT){
						if(bwa_verbose >= 10){
							fprintf(stderr,"Going into timeout mode\n");
							fprintf(stderr,"Starting : %ld\n",qe->starting_read_id);
						}
						vdip = 0xffffffff;
						rc = fpga_pci_poke(fpga_pci_local->pci_bar_handle,0,vdip);
						do { fpga_pci_peek(fpga_pci_local->pci_bar_handle,0,&vled); } while (vled != 0x0);
						time_out = 1;
						break;
					}
				}

				pthread_mutex_unlock (qc->seedex_mut);

				if(time_out == 0){
					f1v.read_right = true;
					read_scores_from_fpga(w, bw_pci_bar_handle,qe,&f1v,0, BATCH_LINE_LIMIT*64*4 + (2*tid+1) * BATCH_LINE_LIMIT/4*64, extension_meta, alnregs);
				}
#else
				read_buffer.clear();
				f1v.read_right = true;
				pthread_mutex_lock (qc->seedex_mut);
				fpga_func_model(w->opt, load_buffer2, load_buffer_entry_idx2, read_buffer);
				pthread_mutex_unlock (qc->seedex_mut);
				// assert(read_buffer.size() == (load_buffer_entry_idx2.size() - 2 + (sizeof(ResultLine::results) / sizeof(ResultEntry))) / (sizeof(ResultLine::results) / sizeof(ResultEntry)) ) ;
				get_all_scores(w,(uint8_t *)read_buffer.data(),read_buffer.size(),qe,&f1v,extension_meta, alnregs);
#endif
			}

			// model validation
			// for(i = 0;i<qe->num;i++){
			// 	for(int j = 0;j<alnregs[i].n;j++){
			// 		mem_alnreg_v * av = &alnregs[i].a[j];
			// 		for (int k = 0; k < av->n; k++) {
			// 			mem_alnreg_t * a = &av->a[k];
			// 			assert(a->score == alnregs_vv[i].a[j].a[k].score);
			// 		}
			// 	}
			// }

			for(i = 0;i<qe->num;i++){
				qe->regs[i] = (mem_alnreg_v *) malloc(sizeof(mem_alnreg_v));
				kv_init(*qe->regs[i]);

				// time_out = 1;
				if(time_out == 1){
					seed_extension(w->opt, w->bwt, w->bns, w->pac, qe->seqs[i]->l_seq, qe->seqs[i]->seq, qe->chains[i], &alnregs[i], &f1, load_buffer1, load_buffer_entry_idx1, load_buffer2, load_buffer_entry_idx2, extension_meta, 0);
					if (alnregs[i].n > 0) {
						kv_copy(mem_alnreg_t, *qe->regs[i], alnregs[i].a[0]);
						kv_destroy(alnregs[i].a[0]);
					}
				} else {
					// Perform postprocess
					for (int j = 0; j < qe->chains[i]->n; ++j) {
						postprocess_alnreg(w->opt, qe->seqs[i]->l_seq, &(qe->chains[i]->a[j]), &(alnregs[i].a[j]), qe->regs[i]);
						kv_destroy(alnregs[i].a[j]);
					}
				}
				mem_alnreg_v * regs = qe->regs[i];
				regs->n = mem_sort_dedup_patch(w->opt, w->bns, w->pac, (uint8_t*)qe->seqs[i]->seq, regs->n, regs->a);

				if (bwa_verbose >= 4) {
					err_printf("* %ld chains remain after removing duplicated chains\n", regs->n);
					for (int ii = 0; ii < regs->n; ++ii) {
						mem_alnreg_t *p = &(regs->a[ii]);
						printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
					}
				}
				for (int ii = 0; ii < regs->n; ++ii) {
					mem_alnreg_t *p = &(regs->a[ii]);
					if (p->rid >= 0 && w->bns->anns[p->rid].is_alt)
						p->is_alt = 1;
				}
				// Free chains now
				free_chains(qe->chains[i]);
				free(alnregs[i].a);
			}

			free(f1v.a);
			free(alnregs);
			f1v.n = 0;

			load_buffer1.clear();
			load_buffer2.clear();
			load_buffer_entry_idx1.clear();
			load_buffer_entry_idx2.clear();
			extension_meta.clear();

		}

		// Grab queue mutex and add queue_element in the queue
		pthread_mutex_lock (q2->mut);
		while (q2->full) {
			if(bwa_verbose >= 18)
				printf ("producer: queue FULL.\n");
			pthread_cond_wait (q2->notFull, q2->mut);
		}
		queueAdd (q2, qe);
		pthread_mutex_unlock (q2->mut);
		pthread_cond_signal (q2->notEmpty);

		if(last_entry){
			break;
		}

	}

	pthread_exit(0);
	//return;
}








void worker2_MT(void *data)
{
	extern int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);
	extern void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a);


	worker2_t *w = (worker2_t *)data;
	queue *q = w->queue1;
	queue_t *qe;

	int last_entry;
	int total_last_entries = 0;
	//worker_t *w = (worker_t*)data;
	int i = 0;


	while(1){
		pthread_mutex_lock (q->mut);
		while (q->empty) {
			if(bwa_verbose >= 18)
				printf (" (T3) queue EMPTY.\n");
			pthread_cond_wait (q->notEmpty, q->mut);
		}
		queueDel (q, &qe);
		pthread_mutex_unlock (q->mut);
		pthread_cond_signal (q->notFull);
		last_entry = 0;
		last_entry = qe->last_entry;
		total_last_entries += last_entry;

		if(last_entry == 0){
			for(i = 0;i<qe->num;i++){
				if (!(w->opt->flag&MEM_F_PE)) {
					if (bwa_verbose >= 4) printf("(T3) =====> Finalizing read '%s' <=====\n", qe->seqs[i]->name);

					mem_mark_primary_se(w->opt, qe->regs[i]->n, qe->regs[i]->a, qe->starting_read_id + i);
					if (w->opt->flag & MEM_F_PRIMARY5) mem_reorder_primary5(w->opt->T, qe->regs[i]);
					mem_reg2sam(w->opt, w->bns, w->pac, qe->seqs[i], qe->regs[i], 0, 0);
					//free(w->regs[i].a);
				} else {
					if (bwa_verbose >= 4) printf("=====> Finalizing read pair '%s' <=====\n", qe->seqs[i<<1|0]->name);
					mem_sam_pe(w->opt, w->bns, w->pac, w->pes, qe->starting_read_id + i, qe->seqs[i<<1], qe->regs[i<<1]);
					//free(w->regs[i<<1|0].a); free(w->regs[i<<1|1].a);
				}
			}
			delete_queue_entry(qe);
		}

		if(total_last_entries >= NUM_FPGA_THREADS){
			delete_queue_entry(qe);
			break;
		}


	}

	pthread_exit(0);
	//return;
}



/*static int
check_slot_config(int slot_id)
{
	int rc;
	struct fpga_mgmt_image_info info = {0};

	// get local image description, contains status, vendor id, and device id
	rc = fpga_mgmt_describe_local_image(slot_id, &info, 0);
	fail_on(rc, out, "Unable to get local image information. Are you running as root?");

	// check to see if the slot is ready
	if (info.status != FPGA_STATUS_LOADED) {
		rc = 1;
		fail_on(rc, out, "Slot %d is not ready", slot_id);
	}

	// confirm that the AFI that we expect is in fact loaded 
	if (info.spec.map[FPGA_APP_PF].vendor_id != pci_vendor_id ||
		info.spec.map[FPGA_APP_PF].device_id != pci_device_id) {
		rc = 1;
		printf("The slot appears loaded, but the pci vendor or device ID doesn't "
			   "match the expected values. You may need to rescan the fpga with \n"
			   "fpga-describe-local-image -S %i -R\n"
			   "Note that rescanning can change which device file in /dev/ a FPGA will map to.\n"
			   "To remove and re-add your edma driver and reset the device file mappings, run\n"
			   "`sudo rmmod edma-drv && sudo insmod <aws-fpga>/sdk/linux_kernel_drivers/edma/edma-drv.ko`\n",
			   slot_id);
		fail_on(rc, out, "The PCI vendor id and device of the loaded image are "
						 "not the expected values.");
	}

out:
	return rc;
}*/


void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0, fpga_pci_conn * fpga_pci)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	extern void kt_for_batch(int n_threads, void (*func)(void*,int,int,int), void *data, int n, int batch_size); // [QA] new kt_for for batch processing
	worker_t w;
	worker2_t w2;
	queue_coll qc[NUM_FPGA_THREADS];
	mem_pestat_t pes[4];
	double ctime, rtime;
	int i,j,total_qual;
	int slot_id = 0;


		ctime = cputime(); rtime = realtime();
		int rc;


	for(i=0;i<n;i++){
		total_qual = 0;
		int l_seq = seqs[i].l_seq;
		for(j=0;j<l_seq;j++){
			total_qual += ((int)seqs[i].qual[j] - 33);
		}
		seqs[i].avg_qual =(uint8_t) (total_qual / l_seq);
	}

		fpga_mem_write_offset = 0;

		fpga_pci_local = fpga_pci;




		global_bns = bns;

		w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
		w.seqs = seqs; w.n_processed = n_processed;
		w.pes = &pes[0];
		w.aux = malloc(opt->n_threads * sizeof(smem_aux_t));
		
		w2.opt = opt; w2.bwt = bwt; w2.bns = bns; w2.pac = pac;
		w2.pes = &pes[0];


		w.n_processed = n;

		// Queue init
		w.queue1 = queueInit();

		if (w.queue1 ==  NULL) {
			fprintf (stderr,"main: Queue Init failed.\n");
			exit (1);
		}

		w2.queue1 = queueInit();
		if (w2.queue1 ==  NULL) {
			fprintf (stderr,"main: Queue 2 Init failed.\n");
			exit (1);
		}

		// SeedEx Mutex
		pthread_mutex_t *seedex_mut = (pthread_mutex_t *) malloc (sizeof (pthread_mutex_t));
		pthread_mutex_init (seedex_mut, NULL);

		for (int j = 0; j < NUM_FPGA_THREADS; ++j) {
			qc[j].q1 = w.queue1;
			qc[j].q2 = w2.queue1;
			qc[j].w = &w;
			qc[j].seedex_mut = seedex_mut;
			qc[j].tid = j;
		}


		for (i = 0; i < opt->n_threads; ++i)
			w.aux[i] = smem_aux_init();
		//kt_for(opt->n_threads, worker1, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
		//kt_for_batch(opt->n_threads, worker1_MT, &w, n, BATCH_SIZE); // find mapping positions
		//
		//

		// Create producer and consumer thread
		

		pthread_t s1, s2[NUM_FPGA_THREADS], s3;
		pthread_create (&s1, NULL, worker1_MT, &w);
		for (int j = 0; j < NUM_FPGA_THREADS; ++j) pthread_create (&s2[j], NULL, fpga_worker, &qc[j]);
		pthread_create (&s3, NULL, worker2_MT, &w2);
		pthread_join (s1, NULL);
		for (int j = 0; j < NUM_FPGA_THREADS; ++j) pthread_join (s2[j], NULL);
		pthread_join (s3, NULL);
		queueDelete (w.queue1);
		queueDelete (w2.queue1);

		pthread_mutex_destroy (seedex_mut);
		free(seedex_mut);


		for (i = 0; i < opt->n_threads; ++i)
			smem_aux_destroy(w.aux[i]);
		free(w.aux);


	/*if (opt->flag&MEM_F_PE) { // infer insert sizes if not provided
		if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
		else mem_pestat(opt, bns->l_pac, n, w.regs, pes); // otherwise, infer the insert size distribution from data
	}*/
	//kt_for(opt->n_threads, worker2, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // generate alignment
	//if (fd >= 0) {
	//close(fd);
	//}
	//rc = fpga_mgmt_close();
	if (bwa_verbose >= 3){
		fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
		fprintf(stderr, " Processed seeds : %ld\n",total_seeds);
	}

}


#ifdef __cplusplus
}
#endif