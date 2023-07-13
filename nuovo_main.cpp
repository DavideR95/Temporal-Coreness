#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <chrono>
#include <utility>
#include <queue>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <signal.h>

#include "temporal_graph.hpp"

// Use (7) if dealing with AS-733 dataset
#define TIME_SLICE (86400 * 7) //(86400 * 5) // (86400*60)

// Uncomment the following define to compile the code for vtune analysis
// #define VTUNE

// Used for vtune code profiling
#if defined(VTUNE) && (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) 
#include "/opt/intel/oneapi/vtune/latest/sdk/include/ittnotify.h"
#endif

// Provare col 90%

/* Likely / unlikely macros for branch prediction */
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

/* Timer 12 hours without vtune, or 15 minutes for vtune */
#ifndef VTUNE
#define TIMEOUT (60 * 60 * 12) // 12 hours
#else
#define TIMEOUT (60 * 1) // 15 mins
#endif

// We assume to work with less than 2^32-1 vertices
using node_t = uint32_t; // This can be changed to uint64_t if needed
using time_type = size_t;

// Timeout handler
volatile bool interrupted = false; 

// Type of a snapshot of a temporal graph
using snapshot_t = std::vector<std::vector<node_t>>; // This is like a static graph

// Heap tree traversal
#define parent(u) ((u-1)/2)
#define left(u) (2*u + 1)
#define right(u) (2*u + 2)

int ceil_log2(unsigned long long x)
{
  static const unsigned long long t[6] = {
    0xFFFFFFFF00000000ull,
    0x00000000FFFF0000ull,
    0x000000000000FF00ull,
    0x00000000000000F0ull,
    0x000000000000000Cull,
    0x0000000000000002ull
  };

  int y = (((x & (x - 1)) == 0) ? 0 : 1);
  int j = 32;
  int i;

  for (i = 0; i < 6; i++) {
    int k = (((x & t[i]) == 0) ? 0 : j);
    y += k;
    x >>= k;
    j >>= 1;
  }

  return y;
}

/* Auxiliary function used to read input graphs */
template <typename T>
temporal_graph_t<T>* ReadGraph(const std::string& input_file, bool is_nde = true) {
    size_t N = 0;
    std::ifstream file(input_file);
    // file >> N;
    // std::vector<std::vector<T>> edg(N);
    // std::vector<std::unordered_map<size_t, std::vector<T>>> temp_edges(N);
    // std::unordered_map<size_t, std::vector<std::vector<T>>> temp_edges;
    std::cerr << "Reading input graph..." << std::endl;
    T u, v;
    time_type timestamp;
    size_t rows = 0;

    std::vector<temporal_edge> edges;
    edges.reserve(10000);
    std::unordered_set<T> vertices;

    if(!file.good()) { std::cerr << "Error opening " << input_file << std::endl; return nullptr; }

    while(!file.eof()) {
        file >> u >> v >> timestamp;
        edges.push_back({u, v, timestamp});
        if(u > N) N = u;
        if(v > N) N = v;
        vertices.insert(u);
        vertices.insert(v);

        // temp_edges[u][timestamp].push_back(v);
        std::cerr << "Righe processate: " << ++rows << "\r";
        //edg[u].push_back(v);
        //edg[v].push_back(u);
    }
    std::cerr << std::endl;


    N++; // !!!!!

    std::cerr << "Max node: " << N << ", actual number of vertices: " << vertices.size() << std::endl;

    // Sort edges based on the timestamp
    std::sort(edges.begin(), edges.end(), compare_edges_timestamp);

    // Group edges according to a time slice (e.g. 1 hour, 1 day, 1 week etc)
    time_type min_timestamp = edges[0].timestamp;
    size_t i = 0;
    for(auto& edge : edges) {
        if(edge.timestamp >= min_timestamp + TIME_SLICE) {
            i++;
            min_timestamp = edge.timestamp;
        }
        edge.timestamp = i; // Rewrite the timestamps starting from 0
    }

    std::cerr << "Max index created: " << i << " (each of " << TIME_SLICE << " seconds)" << std::endl;
    std::cerr << "tau = " << i+1 << std::endl;
    
    return new temporal_graph_t<T>(edges, N, i+1); // new temporal_graph_t<T>(temp_edges, N, max_timestamp - min_timestamp);
}

// Computes the interval of the leaves reachable in a heap
// starting from the internal node tree_node
// depth = max depth of the heap
std::pair<size_t, size_t> leaves(size_t tree_node, size_t depth, size_t tau) {
    int layer = std::log2(tree_node+1);

    // while(true) {
        size_t leaves_start = (tree_node + 1) * (1 << (depth-layer)) - 1; // n * (2^(depth-layer(n))) - 1
        size_t leaves_end = (tree_node + 2) * (1 << (depth-layer)) - 2; // (n+1) * (2^(depth-layer(n))) - 2

        // All of the leaves lie in the second-to-last level of the tree
    //    if(leaves_start > 2 * tau - 2) { 
    //        depth--; // Make it like the tree is one level shorter 
    //        continue; // And redo the computation
    //    }

        return std::make_pair(leaves_start, leaves_end);
    // }

    return std::make_pair(0, 0);
}

// Recursive procedure that merges tree nodes with the same parent
// into that parent. Needed because of the nature of the tree
// that may have leaves in two different layers 
// and the procedure find_tree_nodes is not able to recognize this
std::vector<time_type> fix_interval(std::vector<time_type>& v) {
    if(v.size() <= 1) return v; // Base case, no need to do anything
    bool rec = false; // Should we recur?
    std::vector<time_type> result;
    for(auto i=0;i<v.size()-1;i++) {
        // If two consecutive nodes share the same parent
        if(parent(v[i]) == parent(v[i+1])) {
            result.push_back(parent(v[i])); // Merge them with the parent and put it in the result
            rec = true;
            continue; // Don't put v[i] in the result
        }
        result.push_back(v[i]);
    }

    if(!rec) result.push_back(v.back()); // If we did merge something we may skip the last node, so put it back

    if(rec) return fix_interval(result); // Recur in the new array (go up in the tree)
    else return result;
}

size_t count_parent(size_t node) {
    size_t count = 0;
    while(node > 0) { node = parent(node); count++; }

    return count;
}

std::vector<time_type> find_tree_nodes(temporal_graph_t<node_t>* g, time_type start, time_type end) {



    if(start < 0 || end >= g->tau()) {
        std::cerr << "Invalid time parameters" << std::endl;
        return {};
    }

    std::vector<time_type> result;

    unsigned int v = g->tau(); // compute the next highest power of 2 of 32-bit 
    // Fare un for dove tolgo gli snapshot che sono già empty

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;

    // Align start and end to the actualy nodes in the heap
    start = start + v /* + g->tau()*/ - 1 ;
    end = end + v /* + g->tau() */ - 1;
    size_t a = start, b = end; // start, end will be modified in the loops
    if(a == b) return {start};

    result.resize((b - a) + 1);
    std::iota(result.begin(), result.end(), a);

    return result;

    size_t depth = std::ceil(std::log2(v + g->tau() - 1)) - 1; // Total depth of the tree

    assert(depth > 0);

    // Find out if a and b lie on the same level of the tree
    if(count_parent(a) != count_parent(b)) {
        std::cout << "a = " << a << ", b = " << b << ", tau = " << g->tau() << std::endl;
        assert(false);
        // The tree is unbalanced, we have to break it into two parts
    //     depth--; // Do the upper level first
    //     end = (1 << (int) std::ceil(std::log2(a))) - 2; // The upper level ends at the previous power of 2

    //     size_t curr;
    //     while(start <= end) {
    //         std::pair<size_t, size_t> my_leaves; // Leaves of the subtree rooted in curr
    //         curr = start;
    //         while(true) {
    //             if(curr == 0) break; // We arrived at the root, stop
    //             my_leaves = leaves(parent(curr), depth, g->tau()); // Compute the leaves of curr
    //             if(my_leaves.first < a || my_leaves.second > b) // if [first, last] is not a subset of [a,b]
    //                 break;

    //             curr = parent(curr);   
    //         }            
    //         result.push_back(curr); // current node is part of the solution
    //         // Adjust the interval by doing the set difference between [a,b] \setminus my_leaves
    //         my_leaves = leaves(curr, depth, g->tau());
    //         if(start < my_leaves.first) end = std::min(end, my_leaves.first);
    //         else if(start == my_leaves.first && end >= my_leaves.second) { start = my_leaves.second + 1; }
    //         if(end > my_leaves.second) { start = std::max(start, my_leaves.second); }
    //         else if(end == my_leaves.second && start <= my_leaves.first) { end = my_leaves.first + 1; }
    //     }

    //     // Lower level
    //     depth++;
    //     start = (1 << (int) std::ceil(std::log2(a))) - 1;
    //     a = start;
    //     end = b;
    //     // while(start <= end) {
    //     //     std::pair<size_t, size_t> my_leaves;
    //     //     curr = start;
    //     //     while(true) {
    //     //         if(curr == 0) break;
    //     //         my_leaves = leaves(parent(curr), depth, g->tau());
    //     //         if(my_leaves.first < a || my_leaves.second > b)
    //     //             break;

    //     //         curr = parent(curr);
                
    //     //     }
    //     //     my_leaves = leaves(curr, depth, g->tau());
    //     //     result.push_back(curr);
    //     //     if(start < my_leaves.first) {
    //     //         end = std::min(end, my_leaves.first);
    //     //     }
    //     //     else if(start == my_leaves.first && end >= my_leaves.second) { start = my_leaves.second + 1; }
    //     //     if(end > my_leaves.second) { start = std::max(start, my_leaves.second); }
    //     //     else if(end == my_leaves.second && start <= my_leaves.first) { end = my_leaves.first + 1; }
    //     // }
    // }
    // else {
    //     depth = std::min((double) depth, std::log2(a));
    // }
    }

    size_t curr;
    while(start <= end) {
        std::pair<size_t, size_t> my_leaves;
        curr = start;
        while(true) {
            if(curr == 0) break;
            my_leaves = leaves(parent(curr), depth, g->tau());
            if(my_leaves.first < a || my_leaves.second > b)
                break;

            curr = parent(curr);
            
        }
        my_leaves = leaves(curr, depth, g->tau());
        result.push_back(curr);
        if(start < my_leaves.first) {
            end = std::min(end, my_leaves.first);
        }
        else if(start == my_leaves.first && end >= my_leaves.second) { start = my_leaves.second + 1; }
        if(end > my_leaves.second) { start = std::max(start, my_leaves.second); }
        else if(end == my_leaves.second && start <= my_leaves.first) { end = my_leaves.first + 1; }
    }


    std::sort(result.begin(), result.end());
    return fix_interval(result); // Recursively merge consecutive nodes with the same parents that may have been produced
    // return result;
}

// Intersects two snapshots of a temporal graph
// Done by intersecting the two adjacency lists
snapshot_t snap_intersect(snapshot_t& g, snapshot_t& h) {
    if(g.empty()) return h;
    else if(h.empty()) return g;

    snapshot_t result(g.size());
    for(auto v=0;v<g.size();v++) {
        if(g[v].empty() || h[v].empty()) continue;
        
        std::set_intersection(g[v].begin(), g[v].end(), h[v].begin(), h[v].end(), std::back_inserter(result[v]));
    }

    return result;
}

snapshot_t snap_merge(snapshot_t& g, snapshot_t& h) {
    if(g.empty()) return h;
    else if(h.empty()) return g;
    snapshot_t result(g.size());
    for(auto v=0;v<g.size();v++) {
        assert(std::is_sorted(g[v].begin(), g[v].end()));
        assert(std::is_sorted(h[v].begin(), h[v].end()));
        std::merge(g[v].begin(), g[v].end(), h[v].begin(), h[v].end(), std::back_inserter(result[v]));
        assert(std::is_sorted(result[v].begin(), result[v].end()));
        // size_t i = 0;
        // while(temp[v].size() > 0 && i < temp[v].size()-1) {
        //     size_t count = 1;
        //     while(temp[v][i] == temp[v][i+1]) {
        //         count++;
        //         i++;
        //     }
        //     // Copia temp[v][i] il minimo tra count e min_count.
        //     // Se count < min_count, potrebbe essere utile in una unione successiva.
        //     // Se count >= min_count, non me ne servono di più per farli rimanere nelle prossime unioni
        //     result[v].insert(result[v].end(), count/*std::min(count, min_count)*/, temp[v][i]); 
        //     i++;
        // }
    }

    return result;
}

// Intersect a number of snapshots from the heap tree
// Indexes contains pointers to the nodes in the tree to be intersected
// The result is one static graph made of the intersection of all the tree nodes 
snapshot_t intersect_all(std::vector<snapshot_t>& tree, std::vector<time_type>& indexes) {

    if(indexes.empty()) return {}; // Return an empty snapshot
    else if(indexes.size() < 2) return tree[indexes[0]]; 

    snapshot_t result = tree[indexes[0]];

    for(auto ind=1;ind<indexes.size();ind++) {
        result = snap_intersect(result, tree[indexes[ind]]);
    }

    return result;
}

// Intersect a number of snapshots from the heap tree
// If the last parameter is false, perform the union
// Indexes contains pointers to the nodes in the tree to be intersected
// The result is one static graph made of the intersection/union of all the tree nodes 
snapshot_t union_all(std::vector<snapshot_t>& tree, std::vector<time_type>& indexes, size_t h=1) {

    if(indexes.empty()) return {}; // Return an empty snapshot
    //else if(indexes.size() < 2) return tree[indexes[0]]; 

    snapshot_t result = tree[indexes[0]]; // non aggiungere in indexes gli snapshot che sono vuoti
    snapshot_t final_snapshot(result.size()); // Controllare che non abbia size 0

    // First, merge the snapshots and keep every occurrence of every vertex
    for(auto ind=1;ind<indexes.size();ind++) {
        result = snap_merge(result, tree[indexes[ind]]);
    }

    // Then cleanup the adjacency lists (avoid multiedges)
    for(auto v=0;v<result.size();v++) {
        size_t i = 0;
        while(result[v].size() > 0 && i < result[v].size()) {
            size_t count = 1; // Count each occurrence of a vertex
            while(i+1 < result[v].size() && result[v][i] == result[v][i+1]) {
                count++;
                i++;
            }
            if(count >= h) { // In the final snapshot, put each vertex once (avoid multiedges)
                final_snapshot[v].push_back(result[v][i]);
            }
            i++;
        }
    }

    return final_snapshot;
}

// Build a heap tree off from the timestamps 
std::vector<snapshot_t> build_tree(temporal_graph_t<node_t>* graph, int h = 1) {

    unsigned int v = graph->tau(); // compute the next highest power of 2 of 32-bit v

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;

    std::vector<snapshot_t> tree(v + graph->tau() - 1);

    // std::cout << "tree size: " << tree.size() << std::endl;

    for(auto i=v-1;i<tree.size();i++) {
        tree[i] = graph->get_snapshot(i-v+1); // Timestamps start from zero
        // Leaves in the heap start from index size()/2 
    }

    // Now perform the intersections/union in the parent nodes
    for(int i=v-2;i>=0;i--) { // diventa tau-1 ???
        auto left_child = left(i);
        auto right_child = right(i);
        if(left_child >= tree.size() && right_child >= tree.size()) {
            tree[i] = snapshot_t(); // Quando si fa l'intersezione bisogna accorgersi che questo è vuoto -> ok fatto
        }
        else {
            if(left_child < tree.size() && right_child >= tree.size()) {
                tree[i] = tree[left_child];
            }
            else if(left_child >= tree.size() && right_child < tree.size()) {
                tree[i] = tree[right_child];
            }
            else { // Normal case
                tree[i] = (h < 1) ? snap_intersect(tree[left(i)], tree[right(i)]) : snap_merge(tree[left(i)], tree[right(i)]); // For each internal node
            }
        }   
    }

    // std::vector<snapshot_t> tree(2 * graph->tau() - 1);
    // 2^ceil(log_2(tau)) - 1 è l'indice di partenza
    // E ne devo allocare altri tau -1 dopo Quindi G_tau sta in 2^ceil(log) - 1 + tau -1
    // for(auto i=graph->tau()-1;i<tree.size();i++) { // Leaves
    //     tree[i] = graph->get_snapshot(i-graph->tau()+1); // Timestamps start from zero
    //     // Leaves in the heap start from index size()/2 
    // }

    // // Now perform the intersections/union in the parent nodes
    // for(int i=graph->tau()-2;i>=0;i--) { // diventa tau-1
    //     tree[i] = (h < 1) ? snap_intersect(tree[left(i)], tree[right(i)]) : snap_merge(tree[left(i)], tree[right(i)], h); // For each internal node
    // }

    // std::cout << "L'albero è: " << std::endl;
    // for(auto i=0;i<tree.size();i++) {
    //     if(tree[i].empty()) std::cout << "(" << i << ", __, __) | ";
    //     else { std::cout << "(" << i << ", " << left(i) << ", " << right(i) << ") | "; }
    // }

    return tree;
}

// Finds the k-cores in a snapshot (static graph)
std::vector<node_t> k_cores(snapshot_t& graph, int k) {
    std::vector<int> degree(graph.size());
    node_t max_node = 0;
    size_t max_degree = 0;

    for(auto i=0;i<graph.size();i++) {
        degree[i] = graph[i].size();
        if(degree[i] > max_degree) {
            max_node = i;
            max_degree = degree[i];
        }

    }

    std::queue<node_t> Q; // Will contain the nodes to be removed from the graph
    for(node_t i=0;i<graph.size();i++) {
        if(degree[i] < k)
            Q.push(i); // Push all the nodes with degree less than k
    }

    // While there is a vertex of degree less than k
    while(!Q.empty()) {
        node_t vertex = Q.front();
        Q.pop();
        degree[vertex] = 0;

        // Remove that vertex from the graph
        for(auto& neigh : graph[vertex]) {
            degree[neigh]--;
            if(degree[neigh] > 0 && degree[neigh] < k) {
                Q.push(neigh); // If a neighbor reaches degree less than k, set it up for removal
                degree[neigh] = 0;
            }
        }
    }
    // RIADATTARE PER LA CORENESS DI OGNI NODO
    std::vector<node_t> kcore;
    for(auto v=0;v<graph.size();v++) {
        if(degree[v] >= k) 
            kcore.push_back(v);
    }

    std::cout << "Max degree: " << max_degree << " (node: " << max_node << ")" << std::endl;

    return kcore;
}

std::vector<int> compute_coreness(snapshot_t& graph, size_t& degeneracy) {
    std::vector<std::vector<node_t>> bins(graph.size());
    std::vector<int> degree(graph.size(), 0);
    std::vector<int> coreness(graph.size(), 0);
    size_t max_deg = 0;

    for(auto i=0;i<graph.size();i++) {
        degree[i] = graph[i].size();
        bins[degree[i]].push_back(i);
        coreness[i] = degree[i];
        if(degree[i] > max_deg) max_deg = degree[i];
    }

    size_t k = 1;
    degeneracy = 0;

    while(k <= max_deg) {
        while(k < bins.size() && !bins[k].empty()) {

            node_t v = bins[k].back();
            bins[k].pop_back();

            for(auto& u : graph[v]) {
                if(degree[u] > k) {
                    bins[degree[u]].erase(std::find(bins[degree[u]].begin(), bins[degree[u]].end(), u));
                    degree[u] = (degree[u] == 0) ? 0 : degree[u]-1;
                    if(degree[u] < max_deg && degree[u] > 0)
                        bins[degree[u]].push_back(u);
                }
            }
            coreness[v] = k;
            degeneracy = std::max(k, degeneracy);
        }

        k++;
    }

    return coreness;
}


// Compute the coreness of each vertex in the snapshot passed
std::vector<int> coreness(snapshot_t& graph, size_t& degeneracy) {
    std::vector<int> degree(graph.size());
    
    for(auto i=0;i<graph.size();i++) {
        degree[i] = graph[i].size(); // Initial degree of all nodes
    }

    // Priority queue containing pairs: degree, node
    // pops the minimum element according to the degree
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> Q;

    for(auto i=0;i<graph.size();i++) {
        if(degree[i] != 0)
            Q.push({degree[i], i}); // Push all the vertices in the graph
    }

    std::vector<int> coreness(graph.size(), -1);
    int k = Q.top().first; // This will be the highest k in the graph (the degeneracy)

    while(!Q.empty()) {
        int deg = Q.top().first;
        int v = Q.top().second;
        Q.pop();

        /* if(deg <= k) { // If v (minimum degree) has degree less than the current k
            coreness[v] = k; 
        }
        else {
            k = deg; // Update the degeneracy / maximum degree
            coreness[v] = k;
        }*/
        if(coreness[v] == -1) {
            k = std::max(k, deg);
            coreness[v] = k;

            for(auto u : graph[v]) { // Remove v from the graph
                if(coreness[u] == -1 /* && degree[u] > k */) {
                    degree[u]--;
                    Q.push({degree[u], u});
                }
            }
        }

        // Increment current degree threshold if all vertices with degree >= current threshold have been processed
        // while (!Q.empty() && Q.top().first <= k) {
        //     Q.pop();
        // }
        // k++;

    }
    
    degeneracy = k; // Also return the degeneracy value of the graph (highest k)
    return coreness;
}

int main(int argc, char* argv[]) {
    // Used for vtune code profiling
#if defined(VTUNE) && (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) 
    __itt_pause();
#endif
    bool is_nde = true; // Input format
    unsigned short k = 3;
    int h = 1;
    bool intersect = false;
    time_type window_size;

    if(argc <= 1) {
        // case 3:
        //     k = (unsigned short) std::strtoul(argv[2], NULL, 10);
        //     h = -1;
        //     break;
        // case 4:
        //     k = (unsigned short) std::strtoul(argv[2], NULL, 10);
        //     h = (size_t) std::strtoul(argv[3], NULL, 10);
        //     intersect = false;
        //     break;
        // case 5:
        //     k = (unsigned short) std::strtoul(argv[2], NULL, 10);
        //     h = (size_t) std::strtoul(argv[3], NULL, 10);
        //     intersect = true;
        //     window_size = (time_type) std::strtoul(argv[4], NULL, 10);
        //     break;

        std::cerr << "Usage: " << argv[0] << " <input edge list path> <k> <h = minimum union size>" << std::endl;
        return 1;
    }
    
    // Read the input graph from file
    auto graph = ReadGraph<node_t>(argv[1]);
    if(!graph) {
        std::cerr << "Error while reading the input graph " << argv[1] << ", quitting..." << std::endl;
        return -1;
    }

    // Set the timeout
    signal(SIGALRM, [](int) { interrupted = true; });
    signal(SIGINT, [](int) { interrupted = true; });
    alarm(TIMEOUT);

    // Output stats on the graph (used to understand when the graph is ready)
    std::cerr << "Graph read. " << std::endl;

#if defined(VTUNE) && (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) 
    __itt_resume(); // For vtune profiling
#endif

    bool stop = false;
    int countt = 0;

    // Preprocess the graph to create the trees for both intersection(s) and union(s)
    std::cout << "Building the intersection tree..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto intersection_tree = build_tree(graph, -1);
    auto elapsed_intersection_tree = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Time to build the INT tree: " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_intersection_tree).count()/1000. << " seconds." << std::endl;
    std::cout << "Building the union tree..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto union_tree = build_tree(graph, 1);
    auto elapsed_union_tree = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Time to build the UNION tree: " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_union_tree).count()/1000. << " seconds." << std::endl;
    std::cout << "Done. " << std::endl;

    // Start working
    do {

        std::cout << "Window size? (0 to stop) ";
        std::cin >> window_size;
        if(window_size == 0) { stop = true; break; }
        std::cout << "Union size? (-1 for intersection) ";
        std::cin >> h;
        if(h < 0) { h = -1; intersect = true; }
        else intersect = false;

        if(window_size >= graph->tau()) window_size = graph->tau()-1;
        // Set h to be half the window_size 
        // h = (!intersect) ? window_size / 2 : -1;
        auto start = std::chrono::high_resolution_clock::now();
        
        // auto tree_time = std::chrono::high_resolution_clock::now();

        std::vector<std::vector<size_t>> final_degrees(graph->size());
        std::vector<std::vector<int>> all_coreness(graph->size());

        // std::cout << "H è: " << h << std::endl;

        for(time_type i=0;i+window_size-1<graph->tau();i++) {
            auto tree_nodes_to_use = find_tree_nodes(graph, i, std::min(i+window_size-1, graph->tau()));
            // std::cout << "Snapshots from " << i << " to " << std::min(i+window_size-1, graph->tau()) << std::endl;

            auto find_tree_nodes_time = std::chrono::high_resolution_clock::now();

            snapshot_t final_snapshot;
            
            if(intersect) {
                final_snapshot = intersect_all(intersection_tree, tree_nodes_to_use);
            }
            else {
                final_snapshot = union_all(union_tree, tree_nodes_to_use, h);
            }

            auto intersect_time = std::chrono::high_resolution_clock::now();

            // std::vector<node_t> k_core;

            // if(!final_snapshot.empty()) 
            //     k_core = k_cores(final_snapshot, k);
            auto k_core_time = std::chrono::high_resolution_clock::now();

            size_t degeneracy = 0;
            auto vertex_coreness = compute_coreness(final_snapshot, degeneracy);

            auto total_elapsed = std::chrono::high_resolution_clock::now() - start;

            if(final_snapshot.empty()) { 
                // std::cout << "The resulting graph is empty." << std::endl; countt++; 
                for(auto j=0;j<graph->size();j++) {
                    all_coreness[j].push_back(-1);
                    final_degrees[j].push_back(-1);
                }
            }
            else {
                assert(vertex_coreness.size() == graph->size());
            }

            // std::cout << "Nodi dell'albero intersecati: ";
            // for(auto& i : tree_nodes_to_use) std::cout << i << " ";
            // std::cout << std::endl;

            // for(auto j=0;j<final_snapshot.size();j++) final_degrees[i][j] = final_snapshot[i].size();

            // std::cout << k << "-core: ";
            // for(auto& i : k_core) std::cout << i << " (deg = " << final_snapshot[i].size() << ") ";
            // std::cout << std::endl << std::endl;

            // std::cout << "Degeneracy of the snapshot: " << degeneracy << std::endl << std::endl;

            // std::cout << "Coreness of all nodes: ";
            // fflush(stdout);
            for(auto j=0;j<vertex_coreness.size();j++) {
                // if(vertex_coreness[i] != 0 && final_snapshot[i].size() != 0)
                //if(vertex_coreness[j] > 0 && vertex_coreness[j] != final_snapshot[j].size())
                //if(vertex_coreness[i] != 0)
                //    std::cout << j << " (cness: " << vertex_coreness[j] << ", deg: " << final_snapshot[j].size() << ") ";
                
                all_coreness[j].push_back(vertex_coreness[j]);
                final_degrees[j].push_back(final_snapshot[j].size());
            }
            std::cout << std::endl;

            
            //std::cout << "Time to identify the tree nodes: " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(find_tree_nodes_time-tree_time).count()/1000. << " seconds." << std::endl;
            //std::cout << "Time to intersect/union the snapshots: " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(intersect_time-find_tree_nodes_time).count()/1000. << " seconds." << std::endl;
            // std::cout << "Time to compute the k-core: " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(k_core_time-intersect_time).count()/1000. << " seconds." << std::endl;
            std::cout << "Time to compute the coreness of all nodes: " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(total_elapsed+start-k_core_time).count()/1000. << " seconds." << std::endl;
            std::cout << "Total elapsed time: " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(total_elapsed).count()/1000. << " seconds." << std::endl;

        }

        // std::cout << "Alla fine countttt è " << countt << std::endl;
#if defined(VTUNE) && (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) 
        __itt_detach(); // For vtune profiling
#endif

        // Create a csv with all info
        // #ifdef AAAA
        std::string path(argv[1]);
        std::string base_filename = path.substr(path.find_last_of("/") + 1);
        std::string::size_type const p(base_filename.find_last_of('.'));
        std::string file_without_extension = base_filename.substr(0, p);
        std::cerr << "Il filename è: " << file_without_extension + "_" + std::to_string(graph->tau()) + "_" + std::to_string(window_size) + "_" + std::to_string(h) << std::endl;
        std::ofstream deg_csv("results/" + file_without_extension + "_" + std::to_string(graph->tau()) + "_" + std::to_string(window_size) + "_" + std::to_string(h) + "_deg.csv", std::ios_base::out | std::ios_base::trunc);
        std::ofstream coreness_csv("results/" + file_without_extension + "_" + std::to_string(graph->tau()) + "_" + std::to_string(window_size) + "_" + std::to_string(h) + "_coreness.csv", std::ios_base::out | std::ios_base::trunc);
        
        deg_csv << "node,";
        coreness_csv << "node,";
        for(time_type i=1;i+window_size-1<graph->tau()+1;i++) { 
            deg_csv << "d_" << i << ",";
            coreness_csv << "c_" << i << ",";
        }

        deg_csv.seekp(-1, std::ios_base::cur);
        deg_csv << std::endl;
        coreness_csv.seekp(-1, std::ios_base::cur);
        coreness_csv << std::endl;
        std::cout << "Alla fine i final_degrees hanno dimensione: " << final_degrees.size() << std::endl; 
        for(auto i=0;i<final_degrees.size();i++) {
            deg_csv << i << ",";
            coreness_csv << i << ",";
            if(final_degrees[i].empty()) {
                assert(false);
                for(auto v=0;v<graph->size();v++)
                    deg_csv << 0 << ",";
            }
            else {
                for(auto& deg : final_degrees[i])
                    deg_csv << deg << ",";
            }
            //std::cout << std::endl;
            if(all_coreness[i].empty()) {
                assert(true);
                for(auto v=0;v<graph->size();v++)
                    coreness_csv << 0 << ",";
            }
            else {
                
                for(auto& cness : all_coreness[i])
                    coreness_csv << cness << ",";
            }

            deg_csv.seekp(-1, std::ios_base::cur);
            deg_csv << std::endl;
            coreness_csv.seekp(-1, std::ios_base::cur);
            coreness_csv << std::endl;
        }

        // Matrice degrees[i][j] dove ij = nodo i, grado al passo j
        // Matrice coreness[i][j] dove ij = nodo i, coreness al passo j
        // Stampare come:
        // nodo 0, grado, grado, grado...
        // ....
        // nodo 0, coreness, coreness, coreness... per ogni window_size 

        deg_csv.close();
        coreness_csv.close();
        // #endif
        all_coreness.clear();
        final_degrees.clear();


    } while(!stop);

    delete graph;

    return (interrupted ? 14 : 0); // Exit status will be 14 if the timeout occurred during execution

}