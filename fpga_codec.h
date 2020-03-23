#pragma once

#define PACKED __attribute__((__packed__))

#define BATCH_LINE_LIMIT	16384*2

#define PACKET_MIDDLE       0u
#define PACKET_START        1u
#define PACKET_END          2u
#define PACKET_COMPLETE     3u

#define C_NULL              0x5
#define C_PADDING           0x6

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
    uint8_t lscore;
    uint8_t gscore;
    uint8_t qle;
    uint8_t tle;
    uint8_t gtle;
    uint8_t exception;
    uint8_t spacing[2];
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


#ifndef NO_BWA
#include <array>
#include <stack>
#include <vector>

#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"
#include "bwamem.h"
typedef struct {
        // TODO add alignment entries
        uint8_t score;
        int fpga_entry_present;

} fpga_data_out_t;


template <typename T, size_t N>
struct arraystack
{
    arraystack() : n(0) {};
    void push_back(T&& item) {assert(n < BATCH_LINE_LIMIT); a[n++] = item;}
    T& back() {return (n > 0)? a[n-1] : a[0];}
    T& at(size_t i) {assert(i >= 0 && i < n && "range error"); return a[i];}
    T* begin() {return &a[0];}
    T* end() {return &a[n];}
    size_t size() {return n;}
    T* data() {return a.data();}
    void clear() {n = 0;}
    void reserve(size_t sz) {assert(sz <= N && "size error");}
    std::array<T,N> a;
    size_t n;
};

typedef arraystack<union SeedExLine, BATCH_LINE_LIMIT> LoadBufferTy;
typedef arraystack<union SeedExLine*, BATCH_LINE_LIMIT> LoadBufferPtrTy;
// typedef std::vector<union SeedExLine> LoadBufferTy;
// typedef std::vector<union SeedExLine*> LoadBufferPtrTy;

typedef struct {
    size_t n,m;
    fpga_data_out_t *a;
	bool read_right;
	LoadBufferTy* load_buffer1;
	LoadBufferPtrTy* load_buffer_entry_idx1;
	LoadBufferTy* load_buffer2;
	LoadBufferPtrTy* load_buffer_entry_idx2;
	std::vector<struct extension_meta_t>* extension_meta;
} fpga_data_out_v;
#endif