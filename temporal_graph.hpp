#include <algorithm>
#include <memory>
#include <iostream>
#include <assert.h>
#include <bitset>
#include <iostream>
#include <unordered_map> 

#ifdef __APPLE__
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/avx512.h"
#else
#include <immintrin.h>
#include <malloc.h>
#endif

#include "span.hpp"
#include "cuckoo.hpp"
//#include "./pthash/include/pthash.hpp"

/* Likely / unlikely macros for branch prediction */
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

using node_t = uint32_t;
using time_type = size_t;

#define ALIGNED(x, y) (__builtin_assume_aligned(x, y))

/* Likely / unlikely macros for branch prediction */
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

#include <stdint.h>
static inline uint32_t log2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}

static inline size_t preceding_pow2(size_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return (++v) >> 1;
}

static bool useHashset = false;

inline bool my_binary_search16(node_t* arr, size_t size, const node_t& key) {
    __m512i needle = _mm512_set1_epi32(key);
    __m512i chunk = _mm512_load_si512((__m512i *)arr);
    return _mm512_mask_cmpeq_epu32_mask((1 << size) - 1, chunk, needle) > 0;
}

template <class T>
inline bool my_binary_search32(const T* arr, size_t length, const T& value) {
      long sz = length;
    __m512i needle = _mm512_set1_epi32(value);
    bool result = false;
    __mmask16 mask;
    __mmask16 filter = (1 << 16) - 1;
    while(sz > 0 && !result) {
        __m512i chunk = _mm512_load_si512((__m512i *)arr);
        if(sz < 16) filter = (1 << sz) - 1;//(1 << 16) - 1 - ((1 << (16-sz)) - 1);
        mask = _mm512_mask_cmpeq_epu32_mask(filter, chunk, needle);
        arr += 16;
        sz -= 16;
        result = mask > 0;
    }
    return result;
}

template <class T>
bool my_binary_search(T* arr,/* size_t beg, size_t end,*/size_t length, const T& value) {
    const T* base = arr;
    while (length > 1 && *base != value) { // Branch-less binary partition
      const size_t mid = length/2;
      length -= mid;
      __builtin_prefetch(&base[length / 2 - 1]);
      __builtin_prefetch(&base[mid + length / 2 - 1]);
      base += (base[mid - 1] < value) * mid;
    }
    return *base == value;
}


struct node_struct_t {
    size_t degree;
    node_t* neighs;
    node_t middle;
} __attribute__((aligned(64)));

typedef struct {
    node_t u;
    node_t v;
    time_type timestamp;
} temporal_edge;

bool compare_edges_timestamp(const temporal_edge& a, const temporal_edge& b) {
    return a.timestamp < b.timestamp;
}

template <typename node_t> class temporal_graph_t {
    public:
    temporal_graph_t(std::vector<temporal_edge>& edges, size_t N_, size_t T_, bool is_directed = false) : size_(N_), 
    aligned_edges(N_),
    total_edges_(N_),
    tau_(T_) {

        for(auto& adj : aligned_edges)
            adj.resize(T_);

        for(auto& edge : edges) {
            aligned_edges[edge.u][edge.timestamp].push_back(edge.v);
            // UNDIRECTED GRAPH
            if(!is_directed) aligned_edges[edge.v][edge.timestamp].push_back(edge.u);
        }

        // Sort adj lists in ascending order
        for(auto& v : aligned_edges) {
            for(auto& adj_list : v) {
                std::sort(adj_list.begin(), adj_list.end());
                adj_list.erase(std::unique(adj_list.begin(), adj_list.end()), adj_list.end());
            }
        }

        
    }

    ~temporal_graph_t() {}

    size_t size() const {
        return size_;
    }

    std::vector<std::vector<node_t>> get_snapshot(time_type timestamp) {
        std::vector<std::vector<node_t>> snapshot(size());
        assert(timestamp < tau_);
        for(size_t u=0;u<size();u++) {
            for(auto neighbor : aligned_edges[u][timestamp]) {
                snapshot[u].push_back(neighbor);
            }
        }

        return snapshot;
    }

    time_type tau() const {
        return tau_;
    }

    std::vector<node_t> neighs(node_t v, time_type timestamp = 0) {
        return aligned_edges[v][timestamp];

    }

    nonstd::span<const node_t> neighs_span(node_t v, time_type timestamp = 0) const {
      return nonstd::span<const node_t>(aligned_edges[v][timestamp].data(), aligned_edges[v][timestamp].size());
    }

    //nonstd::span<const node_t> neighs_span(node_t v) const { return nonstd::span<const node_t>(aligned_edges[v].neighs, aligned_edges[v].degree); }
    // node_struct_t neighs(node_t v) const { return aligned_edges[v]; }
    //std::vector<node_t, AlignmentAllocator<node_t, 64>> fwd_neighs(node_t v) const { return edges_[v]; }
    /* nonstd::span<const node_t> fwd_neighs(node_t v, size_t t = 0) const { 
        auto beg = std::upper_bound(aligned_edges[v][t].neighs, aligned_edges[v][t].neighs + aligned_edges[v][t].degree, v);
        auto end = aligned_edges[v][t].degree;
        return nonstd::span<const node_t>(beg, aligned_edges[v][t].neighs + aligned_edges[v][t].degree - beg);//.data(), edges_[v].end());
    } */

    bool are_neighs(node_t a, node_t b, time_type timestamp = 0) {
        return std::find(aligned_edges[a][timestamp].begin(), aligned_edges[a][timestamp].end(), b) != aligned_edges[a][timestamp].end();
    }

    size_t degree(node_t v, time_type timestamp = 0) const {
        return aligned_edges[v][timestamp].size();
    }

    private:
        std::vector<std::vector<std::vector<node_t>>> aligned_edges;
        //std::vector<std::unordered_map<size_t, node_struct_t>> aligned_edges_old;
        time_type tau_; // Number of timestamps
        // std::vector<cuckoo_hash_set<node_t>> hash_edges;
        size_t size_; // Number of vertices
        size_t total_edges_; // Number of all edges in all snapshots
        // size_t avg_list_length, num_are_neighs_calls, num_bin_search, num_less_than_16, avg_length_bin_search, avg_jmp_while;
    
};

