#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <emmintrin.h>

#include <iostream>
#include <verilated.h>          // Defines common routines
#include "Vscore_matrix_compute.h"           // From Verilating "bin2bcd.v"
#include "sys_defs.h"

#define C_NULL 5
#define C_FILL 6  // filler string to drain buffer

#ifdef __cplusplus
extern "C" {
#endif

int sequence_counter = 0;

int get_w (int qlen, int w, const int8_t *mat)
{
    int i, k, max, max_ins, max_del;
    for (i = 0, max = 0; i < 25; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + 5 - 6) / 1 + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + 5 - 6) / 1 + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del; // TODO: is this necessary?
    return w;
}

int ksw_extend_hw(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
    Vscore_matrix_compute *top = new Vscore_matrix_compute();       // Create local instance for multi threading

    const uint8_t *query_itr, *target_itr;
    int score = 0;
    int row_score = 0;
    int gscore = 0;
    int row_gscore = 0;
    int max_i = 0, max_j = 0, max_ie = 0;
    int _max_i = 0, _max_j = 0, _max_ie = 0;
    int curr_i = 1;

    unsigned int main_time = 0;     // Current simulation time
    bool computing = false;
    bool compute_done = false;
    int min_seq_len = (int)(sizeof(top->score_matrix_compute__DOT__h_row1) * 2 + 2) ; // BW*2 to drain...
    int num_chars_sent = 0;
    bool zero_row_exists = false, row_score_released = false;
    unsigned int start_time = 0;
    int initial_g_skip = 0;
    int remaining_gscore_fetch = 0;
    bool extension_exception = false;

    top->rst = 1;                         // Set some inputs
    top->clk = 0;
    top->start_i = 0;
    top->init_score_i = h0;
    top->query_entry_i = C_NULL;
    top->target_entry_i = C_NULL;


    while (!Verilated::gotFinish()) {

        if ((main_time % 5) == 0)
            top->clk = !top->clk;         // Toggle clock

        if (main_time > 10) {             // Release reset
            top->rst = 0;
        }
                                          // Assert start flag
        if (main_time > 21 && !computing && top->start_i == 0) {
            top->start_i = 1;
            // printf("#%d setting start_i\n", main_time);
            start_time = main_time;
        }


                                          // Query sequence
        if ((main_time % 10) == 1 && (top->start_i == 1 || computing)) {
            if (!computing) {
                // printf("#%d sequencing start\n", main_time);
                top->start_i = 0;
                computing = true;
                query_itr = query;
                target_itr = target;

                union {unsigned int u[4]; unsigned long l[2];} pack1, pack2;
                w = 1 + get_w(qlen, w, mat);
                pack1.l[0] = (BW/2 + (w+1)/2 >= 64)? -1 : ((1ul << (BW/2 + (w+1)/2)) - 1);
                pack1.l[1] = (BW/2 + (w+1)/2 >= 64)? ((1ul << (BW/2 + (w+1)/2) - 64) - 1) : 0;
                pack2.l[0] = (BW/2 + (w)/2 >= 64)? -1 : ((1ul << (BW/2 + (w)/2)) - 1);
                pack2.l[1] = (BW/2 + (w)/2 >= 64)? ((1ul << (BW/2 + (w)/2) - 64) - 1) : 0;
                for (int i = 0; i < 4; ++i) top->pe1_en_i[i] = pack1.u[i];
                for (int i = 0; i < 4; ++i) top->pe2_en_i[i] = pack2.u[i];
            }

            top->query_entry_i = C_NULL;
            top->target_entry_i = C_NULL;
            if (query_itr - query < qlen || target_itr - target < tlen) {
                top->query_entry_i = C_FILL;
                top->target_entry_i = C_FILL;
                if (query_itr - query < qlen) {
                    top->query_entry_i = *query_itr;
                    query_itr++;
                }
                if (target_itr - target < tlen) {
                    top->target_entry_i = *target_itr;
                    target_itr++;
                }
                num_chars_sent++;
                // printf("#%d sequencing char (%c,%c) at %d\n", main_time, "ACGTN"[(int)*(query_itr-1)], "ACGTN"[(int)*(target_itr-1)], query_itr - query);
            } else if (num_chars_sent < min_seq_len) {
                top->query_entry_i = C_FILL;
                top->target_entry_i = C_FILL;
                num_chars_sent++;
            }
        }

        top->eval();                      // Evaluate model

                                            // Monitor exception
        if ((main_time % 10) == 1 && computing && top->extension_exception_o) {        // Monitor exception
            extension_exception = true;
        }
                                            // Monitor score output
        if ((main_time % 10) == 1 && computing && (main_time - start_time)/10 > BW) {
            if (!row_score_released) {
                // initial_g_skip = (qlen <= BW) ? BW/2 - qlen/2 + 1 : 0;
                initial_g_skip = 0;
                remaining_gscore_fetch = tlen < qlen * 2 ? tlen : qlen * 2;
                row_score = top->init_score_i;
            }
            max_i = row_score < top->max_row_score_o? curr_i : max_i;
            max_j = row_score < top->max_row_score_o? top->max_row_j_idx_o : max_j;
            row_score = row_score < top->max_row_score_o? top->max_row_score_o: row_score;
            if (initial_g_skip > 0) {
                initial_g_skip--;
            } else if (remaining_gscore_fetch > 0) {
                int max_row_gscore;
                max_row_gscore = top->row_gscore_o;
                remaining_gscore_fetch--;

                max_ie = row_gscore <= max_row_gscore ? curr_i : max_ie;
                row_gscore = row_gscore < max_row_gscore ? max_row_gscore : row_gscore;
            }

            if (top->max_row_score_o > 0 && ! row_score_released) {
                row_score_released = true;
            } else if (top->max_row_score_o == 0 && row_score_released && !zero_row_exists) {
                zero_row_exists = true;
                score = row_score < h0 ? h0 : row_score;
                gscore = row_gscore;
                _max_i = row_score < h0 ? curr_i : max_i;
                _max_j = row_score < h0 ? top->max_row_j_idx_o : max_j;
                _max_ie = max_ie;
            }

            curr_i++;
        }

        //                                     // Monitor score output
        // if ((main_time % 10) == 1 && computing && top->max_row_score_o >= 0) {
        //     // row_score = row_score < top->max_score_o? top->max_score_o: row_score;
        //     row_score = row_score < top->max_row_score_o? top->max_row_score_o: row_score;
        //     row_gscore = row_gscore < top->max_gscore_o? top->max_gscore_o: row_gscore;
        //     if (top->max_row_score_o > 0 && ! row_score_released) {
        //         row_score_released = true;
        //     } else if (top->max_row_score_o == 0 && row_score_released && !zero_row_exists) {
        //         zero_row_exists = true;
        //         score = row_score < h0 ? h0 : row_score;
        //         gscore = row_gscore;
        //     }
        // }

                                          // Wait for done
        if ((main_time % 10) == 1 && computing && top->done_o) {
            compute_done = true;
        }

                                          // Read the output
        else if ((main_time % 10) == 1 && compute_done) {
            // score = top->max_score_o;
            score = !score? row_score : score;
            gscore = !zero_row_exists? row_gscore : gscore;
            max_i = !zero_row_exists? max_i : _max_i;
            max_j = !zero_row_exists? max_j : _max_j;
            max_ie = !zero_row_exists? max_ie : _max_ie;
            // printf("#%d\t SEQ %d\t initial_score %d\t computed_score %d computed_gscore: %d\n", main_time, sequence_counter, h0, score, gscore);
            computing = false;
            compute_done = false;
            num_chars_sent = 0;
            zero_row_exists = false;
            row_score_released = false;
            sequence_counter++;
            if (extension_exception) {
                gscore = gscore * -1 + -1000;
            }
            extension_exception = false;
            // printf("@@@ Done processing\n");
            break;
        }

        main_time++;                      // Time passes...

    }

    *_gscore = gscore;
    *_qle = max_j;
    *_tle = max_i;
    *_gtle = max_ie;

                                          // Read the output
    top->final();

    delete top;

    return score;
}

#ifdef __cplusplus
}
#endif
