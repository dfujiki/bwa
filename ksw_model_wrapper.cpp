#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <emmintrin.h>

#include <iostream>
#include <verilated.h>          // Defines common routines
#include "Vscore_matrix_compute.h"           // From Verilating "bin2bcd.v"

#define C_NULL 5
#define C_FILL 6  // filler string to drain buffer

#ifdef __cplusplus
extern "C" {
#endif

int sequence_counter = 0;

int ksw_extend_hw(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
    Vscore_matrix_compute *top = new Vscore_matrix_compute();       // Create local instance for multi threading

    const uint8_t *query_itr, *target_itr;
    int score = 0;
    int row_score = 0;

    unsigned int main_time = 0;     // Current simulation time
    bool computing = false;
    bool compute_done = false;
    int min_seq_len = (int)(sizeof(top->score_matrix_compute__DOT__h_row1) * 2) ; // BW*2 to drain...
    int num_chars_sent = 0;
    bool zero_row_exists = false, row_score_released = false;

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
        }


                                          // Query sequence
        if ((main_time % 10) == 1 && (top->start_i == 1 || computing)) {
            if (!computing) {
                // printf("#%d sequencing start\n", main_time);
                top->start_i = 0;
                computing = true;
                query_itr = query;
                target_itr = target;
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

                                            // Monitor score output
        if ((main_time % 10) == 1 && computing && top->max_row_score_o >= 0 && !zero_row_exists) {
            // row_score = row_score < top->max_score_o? top->max_score_o: row_score;
            row_score = row_score < top->max_row_score_o? top->max_row_score_o: row_score;
            if (top->max_row_score_o > 0 && ! row_score_released) {
                row_score_released = true;
            } else if (top->max_row_score_o == 0 && row_score_released) {
                zero_row_exists = true;
                score = row_score < h0 ? h0 : row_score;
            }
        }

                                          // Wait for done
        if ((main_time % 10) == 1 && computing && top->done_o) {
            compute_done = true;
        }

                                          // Read the output
        else if ((main_time % 10) == 1 && compute_done) {
            // score = top->max_score_o;
             score = !score? row_score : score;
            // printf("#%d\t SEQ %d\t initial_score %d\t computed_score %d\n", main_time, sequence_counter, h0, top->max_score_o);
            computing = false;
            compute_done = false;
            num_chars_sent = 0;
            zero_row_exists = false;
            row_score_released = false;
            sequence_counter++;
            // printf("@@@ Done processing\n");
            break;
        }

        main_time++;                      // Time passes...

    }
                                          // Read the output
    top->final();

    delete top;

    return score;
}

#ifdef __cplusplus
}
#endif
