#ifndef _graphgenBP_h
#define _graphgenBP_h
#include <vector>

typedef struct /* the bipartite graph data structure */
{
    long n;            // numver of vertices in both sides
    long nrows;        // number of vertices in the left side
    long m;            // number of edges
    long *vtx_pointer; // an array of size n+1 storing the pointer in endV array  ptrs
    long *endV;        // an array of size m that stores the second vertex of an edge. ids
    double *weight;    // not used in unweighted graph
    long *degree;
} graph;

void process_mtx_compressed(char *fname, graph *bGraph);
void fast_mtx_read_build(const char *fname, graph *bGraph);
bool isEqual(graph *bGraph1, graph *bGraph2);

// 多源bfs
long *MS_BFS_Graft(graph *G, long *mateI);

long *Pothen_Fan_Fairnes(graph *G, long *mateI);

/** Clean up. */
void free_graph(graph *bGraph);
graph *swap_side(graph *bGraph);

#endif
