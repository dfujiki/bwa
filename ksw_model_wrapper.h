
#ifndef __AC_KSW_W_H
#define __AC_KSW_W_H

#include <stdint.h>
int ksw_extend_hw(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);

#endif