#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
#include "memsort.h"

#include <fpga_pci.h>
#include <fpga_mgmt.h>
#include <utils/lcd.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <poll.h>
#include <sys/time.h>

#define	MEM_16G		(1ULL << 34)
#define NUM_BINS    200
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

int fd1 = -1;




typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

static uint16_t pci_vendor_id = 0x1D0F; /* Amazon PCI Vendor ID */
static uint16_t pci_device_id = 0xF000;
uint64_t fpga_mem_write_offset1 = 0;




static int
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
}

int fastmap_write_to_fpga(uint32_t *buffer, int channel, size_t buffer_size){

    int rc = 0;
    rc = pwrite(fd1,
            buffer + 0,
            (buffer_size),
            channel*MEM_16G + fpga_mem_write_offset1
            );
    if (rc < 0) {
        fprintf(stderr,"[ERROR]call to pwrite failed.\n");
    }
    fpga_mem_write_offset1 += rc;
    fprintf(stderr,"rc : %d\n",rc);
    return rc;

}

void fastmap_push_to_lsb_64(uint32_t *ar, int size, uint64_t new_char, int shift){
    int carry = 0;
    int shift_internal = shift;
    int i;
    while (shift_internal--) {                           // For each bit to shift ...
        for (i = 0; i < size; i++) {   // For each element of the array from high to low ...
            int next = (ar[i] & 0x80000000) ? 1 : 0;  // ... if the low bit is set, set the carry bit.
            ar[i] = carry | (ar[i] << 1);       // Shift the element one bit left and addthe old carry.
            carry = next;                       // Remember the old carry for next time.
        }   
    }
    
    uint32_t char1= new_char & 0xFFFFFFFF;
    uint32_t char2= new_char>>32;


    ar[0] = char1 | ar[0];
    ar[1] = char2 | ar[1];
}

    void encode_ann_seqs_data(uint64_t seq_start, uint64_t seq_end, uint64_t seq2_start, uint64_t seq2_end, uint32_t *ar){
        int ar_size = 8;

        // Seq2_end 
        fastmap_push_to_lsb_64(ar,ar_size,seq2_end,64);
        
        // Seq2_start 
        fastmap_push_to_lsb_64(ar,ar_size,seq2_start,64);
        
        // Seq_end 
        fastmap_push_to_lsb_64(ar,ar_size,seq_end,64);
        
        // Seq_start 
        fastmap_push_to_lsb_64(ar,ar_size,seq_start,64);
    }

int write_ann_data_to_fpga(const bntseq_t *bns, int channel){
    int64_t l_pac_2 = (bns->l_pac);
    if(bwa_verbose >= 10){
        printf("l_pac for ann_data_to_fpga : %lu\n",l_pac_2);
    }
    int n_seqs = bns->n_seqs;
    int i = 0;
    int rc = 0;
    size_t ann_write_buffer_size = n_seqs * 32;
    uint32_t *ann_write_buffer = malloc(ann_write_buffer_size); // TODO: Magic number size of each seq entry
    int ann_write_buffer_index = 0;
    for(i=0;i<n_seqs;i++){
        bntann1_t ann = bns->anns[i];
        uint64_t seq_start = ann.offset;
        uint64_t seq_end = ann.offset + ann.len - 1;
        uint64_t seq2_start = l_pac_2 - (ann.offset + ann.len); 
        uint64_t seq2_end = l_pac_2 - ann.offset - 1; 
        if(bwa_verbose >= 10){
            printf("@%lu %lu %lu %lu\n",seq_start,seq_end,seq2_start,seq2_end);
        }
        encode_ann_seqs_data(seq_start, seq_end, seq2_start,seq2_end,ann_write_buffer + ann_write_buffer_index);
        ann_write_buffer_index += 8;
    }

    if(bwa_verbose >= 10){
        printf("Ann write buffer size : %zd\n",ann_write_buffer_size);
    }
    rc = fastmap_write_to_fpga(ann_write_buffer,channel,ann_write_buffer_size);

    if(ann_write_buffer)
        free(ann_write_buffer);

    return rc;
}


int write_pac_to_fpga(const uint8_t *pac, int64_t l_pac,int channel){
    //int64_t l_pac_2 = (l_pac)>>2;
    size_t l_pac_bytes = (l_pac)>>2;
    //printf("For pac writing : %lu\n",l_pac_2);
    //write_buffer_capacity = 256 * 16 * sizeof(uint32_t);

    int rc = fastmap_write_to_fpga((uint32_t*)(pac),channel,l_pac_bytes);
    return 0; 
}


static void *process(void *shared, int step, void *_data)
{
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
		return ret;
	} else if (step == 1) {
		const mem_opt_t *opt = aux->opt;
		//const bwaidx_t *idx = aux->idx;
		bwaidx_t *idx = aux->idx;
		if (opt->flag & MEM_F_SMARTPE) {
			bseq1_t *sep[2];
			int n_sep[2];
			mem_opt_t tmp_opt = *opt;
			bseq_classify(data->n_seqs, data->seqs, n_sep, sep);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] %d single-end sequences; %d paired-end sequences\n", __func__, n_sep[0], n_sep[1]);
			if (n_sep[0]) {
				tmp_opt.flag &= ~MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, n_sep[0], sep[0], 0, fd1);
				for (i = 0; i < n_sep[0]; ++i)
					data->seqs[sep[0][i].id].sam = sep[0][i].sam;
			}
			if (n_sep[1]) {
				tmp_opt.flag |= MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0, fd1);
				for (i = 0; i < n_sep[1]; ++i)
					data->seqs[sep[1][i].id].sam = sep[1][i].sam;
			}
			free(sep[0]); free(sep[1]);
		} else mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seqs, aux->pes0, fd1);
		aux->n_processed += data->n_seqs;
		return data;
	} else if (step == 2) {

		for (i = 0; i < data->n_seqs; ++i) {
            
			//if (data->seqs[i].sam) err_fputs(data->seqs[i].sam, stdout);
            process_sam_record(&data->seqs[i]);
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual); free(data->seqs[i].sam);
		}
		free(data->seqs); free(data);
		return 0;
	}
	return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

int main_mem(int argc, char *argv[])
{
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, ignore_alt = 0, no_mt_io = 0;
	int fixed_chunk_size = -1;
	gzFile fp, fp2 = 0;
	char *p, *rg_line = 0, *hdr_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	mem_pestat_t pes[4];
	ktp_aux_t aux;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	while ((c = getopt(argc, argv, "51qpaMCSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') no_mt_io = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
		else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'u') opt->flag |= MEM_F_XB;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'o' || c == 'f') xreopen(optarg, "wb", stdout);
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'h') {
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		} else if (c == 'I') { // specify the insert size distribution
			aux.pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else return 1;
	}

	if (rg_line) {
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 1 >= argc || optind + 3 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
//		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "\nScoring options:\n\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overridden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
		fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
		fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
		fprintf(stderr, "       -o FILE       sam file to output results to [stdout]\n");
		fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
		fprintf(stderr, "       -5            for split alignment, take the alignment with the smallest coordinate as primary\n");
		fprintf(stderr, "       -q            don't modify mapQ of supplementary alignments\n");
		fprintf(stderr, "       -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
			if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
			if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0) {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1; // FIXME memory leak
		}
	} else update_a(opt, &opt0);
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	aux.idx = bwa_idx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);
	if (ignore_alt)
		for (i = 0; i < aux.idx->bns->n_seqs; ++i)
			aux.idx->bns->anns[i].is_alt = 0;

	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	aux.ks = kseq_init(fp);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 2], &fd2);
			if (ko2 == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
			aux.ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}

    
    int rc = 0;

    rc = fpga_mgmt_init();
    //char *afi_id = "agfi-04679145e6d52c9a0";
    int slot_id = 0;

    char device_file_name[256];
    rc = sprintf(device_file_name, "/dev/edma%i_queue_1", slot_id);
    if(rc<0){
        fprintf(stderr,"Unable to formate device file name");
    }
    //fail_on((rc = (rc < 0)? 1:0), out, "Unable to format device file name.");
    //printf("device_file_name=%s\n", device_file_name);

    // make sure the AFI is loaded and ready
    rc = check_slot_config(slot_id);
    //fail_on(rc, out, "slot config is not correct");

    fd1 = open(device_file_name, O_RDWR);
    if(fd1<0){
        fprintf(stderr,"Cannot open device file %s.\nMaybe the EDMA "
               "driver isn't installed, isn't modified to attach to the PCI ID of "
               "your CL, or you're using a device file that doesn't exist?\n"
               "See the edma_install manual at <aws-fpga>/sdk/linux_kernel_drivers/edma/edma_install.md\n"
               "Remember that rescanning your FPGA can change the device file.\n"
               "To remove and re-add your edma driver and reset the device file mappings, run\n"
               "`sudo rmmod edma-drv && sudo insmod <aws-fpga>/sdk/linux_kernel_drivers/edma/edma-drv.ko`\n",
               device_file_name);
        fprintf(stderr,"Unable to open DMA queue");
        //fail_on((rc = (fd < 0)? 1:0), out, "unable to open DMA queue. ");
    }

    int channel = 2;
    fpga_mem_write_offset1 = 0;
    /*write_ann_data_to_fpga(aux.idx->bns,channel);
    sleep(2);
    fprintf(stderr,"Lseek after ann write : %zd\n",lseek(fd1,0,SEEK_CUR));

    rc = fsync(fd1);
    fprintf(stderr,"Lseek after fsync: %zd\n",lseek(fd1,0,SEEK_CUR));
    if(rc < 0){
        fprintf(stderr,"[ERROR] call to fsync failed");
    }
    write_pac_to_fpga(aux.idx->pac,aux.idx->bns->l_pac,channel);
    fprintf(stderr,"Lseek after write pac: %zd\n",lseek(fd1,0,SEEK_CUR));
    rc = fsync(fd1);
    fprintf(stderr,"Lseek after fsync: %zd\n",lseek(fd1,0,SEEK_CUR));*/
    int64_t l_pac_2 = aux.idx->bns->l_pac;
    aux.idx->bns->l_pac = (l_pac_2>>1);
    sleep(2);

    rc = fsync(fd1);
    if(rc < 0){
        fprintf(stderr,"[ERROR] call to fsync failed");
    }

    // Create bins
    //
    sorting_init(aux.idx->bns->l_pac);

	bwa_print_sam_hdr(aux.idx->bns, hdr_line);
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	//aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size;
	kt_pipeline(no_mt_io? 1 : 2, process, &aux, 3);
	free(hdr_line);
	free(opt);
	bwa_idx_destroy(aux.idx);
	kseq_destroy(aux.ks);
	err_gzclose(fp); kclose(ko);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		err_gzclose(fp2); kclose(ko2);
	}

    
    rc = fpga_mgmt_close();

    if (fd1 >= 0) {
        close(fd1);
    }
    

    struct timeval sort_st, sort_et;
    gettimeofday(&sort_st,NULL);
    get_sorted_sam();
    gettimeofday(&sort_et,NULL);
    double sorting_time = ((double)sort_et.tv_sec - (double)sort_st.tv_sec) + (double)((double)(sort_et.tv_usec - sort_st.tv_usec) / (double)(1000000));
    fprintf(stderr,"Time to create sorted sam file : %f\n",sorting_time);
    sorting_close();



	return 0;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0, min_intv = 1, max_len = INT_MAX;
	uint64_t max_intv = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	smem_i *itr;
	const bwtintv_v *a;
	bwaidx_t *idx;

	while ((c = getopt(argc, argv, "w:l:pi:I:L:")) >= 0) {
		switch (c) {
			case 'p': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'i': min_intv = atoi(optarg); break;
			case 'I': max_intv = atol(optarg); break;
			case 'L': max_len  = atoi(optarg); break;
		    default: return 1;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa fastmap [options] <idxbase> <in.fq>\n\n");
		fprintf(stderr, "Options: -l INT    min SMEM length to output [%d]\n", min_len);
		fprintf(stderr, "         -w INT    max interval size to find coordiantes [%d]\n", min_iwidth);
		fprintf(stderr, "         -i INT    min SMEM interval size [%d]\n", min_intv);
		fprintf(stderr, "         -L INT    max MEM length [%d]\n", max_len);
		fprintf(stderr, "         -I INT    stop if MEM is longer than -l with a size less than INT [%ld]\n", (long)max_intv);
		fprintf(stderr, "\n");
		return 1;
	}

	fp = xzopen(argv[optind + 1], "r");
	seq = kseq_init(fp);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
	itr = smem_itr_init(idx->bwt);
	smem_config(itr, min_intv, max_len, max_intv);
	while (kseq_read(seq) >= 0) {
		err_printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (print_seq) {
			err_putchar('\t');
			err_puts(seq->seq.s);
		} else err_putchar('\n');
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
		while ((a = smem_next(itr)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				if ((uint32_t)p->info - (p->info>>32) < min_len) continue;
				err_printf("EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
				if (p->x[2] <= min_iwidth) {
					for (k = 0; k < p->x[2]; ++k) {
						bwtint_t pos;
						int len, is_rev, ref_id;
						len  = (uint32_t)p->info - (p->info>>32);
						pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(idx->bns, pos, len, &ref_id);
						err_printf("\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
					}
				} else err_puts("\t*");
				err_putchar('\n');
			}
		}
		err_puts("//");
	}

	smem_itr_destroy(itr);
	bwa_idx_destroy(idx);
	kseq_destroy(seq);
	err_gzclose(fp);
	return 0;
}

