#pragma once

#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"
#include "bwamem.h"

#include <vector>

#define PACKED __attribute__((__packed__))

struct PACKED LineParams
{
    uint32_t seq_id;
    uint8_t qlen;
    uint8_t tlen;
    uint8_t init_score;
    uint8_t w;
};


struct PACKED SeedExLineTy1
{
    uint8_t preamble;
    struct LineParams params;
    uint8_t spacing;
    uint8_t payload1[27];
    uint8_t payload2[27];
};

struct PACKED SeedExLineTy0
{
    uint8_t preamble;
    uint8_t payload[63];
};

struct PACKED ResultEntry
{
    uint32_t seq_id;
    int8_t lscore;
    int8_t gscore;
    uint8_t qle;
    uint8_t tle;
    uint8_t gtle;
    uint8_t exception;
    uint8_t spacing;
};

struct PACKED ResultLine
{
    uint8_t preamble[4];
    struct ResultEntry results[5];
};

union SeedExLine
{
    struct SeedExLineTy1 ty1;
    struct SeedExLineTy0 ty0;
	struct ResultLine ty_r;
};

struct extension_meta_t
{
	uint32_t read_idx;
	uint32_t chain_id;
	uint32_t seed_id;
};


typedef struct {
        // TODO add alignment entries
        uint8_t score;
        int fpga_entry_present;

} fpga_data_out_t;


typedef struct {
    size_t n,m;
    fpga_data_out_t *a;
	bool read_right;
	std::vector<union SeedExLine>* load_buffer1;
	std::vector<union SeedExLine*>* load_buffer_entry_idx1;
	std::vector<union SeedExLine>* load_buffer2;
	std::vector<union SeedExLine*>* load_buffer_entry_idx2;
	std::vector<struct extension_meta_t>* extension_meta;
} fpga_data_out_v;
