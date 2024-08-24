#include "../include/dgmvm.h"
#include "../include/kasi.h"

// DFS with lookahead that finds a single augmenting path
// called from pothen-fan
long findAugPathLookahead(long sFirst, long *flagLookahead, long *flagDFS, long *mate, graph *G, long *path, long *edgeIndexLookahead, long *edgeIndex)
{

    long *edgeStart = G->vtx_pointer;
    long *endVertex = G->endV;
    long NE = G->m;
    long top = -1;
    path[++top] = sFirst; // push , path is equivalent to stack

    while (top >= 0) // while stack not empty
    {
        long u = path[top];
        long uDegree = edgeStart[u + 1] - edgeStart[u];
        // lookahed part
        while (++edgeIndexLookahead[u] < uDegree)
        {
            long v = endVertex[edgeStart[u] + edgeIndexLookahead[u]];
            if (__sync_fetch_and_add(&flagLookahead[v], 1) == 0)
            {
                if (mate[v] == -1)
                {
                    __sync_fetch_and_add(&flagDFS[v], 1);
                    path[++top] = v; // push
                    return top + 1;  // top = augmenting path length
                }
            }
        }

        while (++edgeIndex[u] < uDegree)
        {

            long v = endVertex[edgeStart[u] + edgeIndex[u]];
            if (__sync_fetch_and_add(&flagDFS[v], 1) == 0)
            {
                if (mate[v] != -1) // means other vertex already allocate this in lookahed phase
                {

                    path[++top] = v;       // push v
                    path[++top] = mate[v]; // push next u
                    break;
                }
            }
        }
        if (edgeIndex[u] == uDegree)
        {
            top -= 2; // pop
        }
    }
    return top + 1;
}

// DFS with lookahead that finds a single augmenting path
// called from pothen-fan
long findAugPathLookaheadReverse(long sFirst, long *flagLookahead, long *flagDFS, long *mate, graph *G, long *path, long *edgeIndexLookahead, long *edgeIndex)
{

    long *edgeStart = G->vtx_pointer;
    long *endVertex = G->endV;
    long NE = G->m;
    long top = -1;
    path[++top] = sFirst; // push , path is equivalent to stack

    while (top >= 0) // while stack not empty
    {
        long u = path[top];
        long uDegree = edgeStart[u + 1] - edgeStart[u];
        // lookahed part
        while (++edgeIndexLookahead[u] < uDegree)
        {
            long v = endVertex[edgeStart[u] + edgeIndexLookahead[u]];
            if (__sync_fetch_and_add(&flagLookahead[v], 1) == 0)
            {
                if (mate[v] == -1)
                {
                    __sync_fetch_and_add(&flagDFS[v], 1);
                    path[++top] = v; // push
                    return top + 1;  // top = augmenting path length
                }
            }
        }

        // while(++edgeIndex[u] < uDegree)
        while (--edgeIndex[u] >= 0)
        {

            long v = endVertex[edgeStart[u] + edgeIndex[u]];
            if (__sync_fetch_and_add(&flagDFS[v], 1) == 0)
            {
                if (mate[v] != -1) // means other vertex already allocate this in lookahed phase
                {

                    path[++top] = v;       // push v
                    path[++top] = mate[v]; // push next u
                    break;
                }
            }
        }
        if (edgeIndex[u] == -1)
        {
            top -= 2; // pop
        }
    }
    return top + 1;
}

// ------------- PF with Fairness ---------------------

long *Pothen_Fan_Fairnes(graph *G, long *mateI)
{

    double time2, time;
    // time = omp_get_wtime();
    long NE = G->m;
    long NV = G->n;
    long *endVertex = G->endV;
    long *edgeStart = G->vtx_pointer;

    long *unmatchedU = (long *)malloc(NV * sizeof(long));
    long *tQ = (long *)malloc(NV * sizeof(long));
    long *mate = (long *)malloc(NV * sizeof(long));
    long *flagDFS = (long *)malloc(NV * sizeof(long));
    long *flagLookahead = (long *)malloc(NV * sizeof(long));
    long *edgeIndex = (long *)malloc(NV * sizeof(long));
    long *edgeIndexLookahead = (long *)malloc(NV * sizeof(long));

    // 记录未匹配的顶点
    long numUnmatchedU = 0;
    for (long u = 0; u < G->nrows; u++)
    {
        if (mateI[u] == -1 && (edgeStart[u + 1] > edgeStart[u]))
        {
            unmatchedU[numUnmatchedU++] = u;
        }
    }

    for (long i = 0; i < NV; i++)
    {
        mate[i] = mateI[i];
        flagLookahead[i] = 0;
        edgeIndexLookahead[i] = -1;
    }

    // 记录增广路径
    long *augPath = (long *)malloc(NV * sizeof(long));

    long iterations = 0;
    while (1)
    {
        iterations++;
        long tQ_len = 0;
        if (iterations % 2 == 1) // 奇数轮更新相关
        {
            for (long i = 0; i < NV; i++)
            {
                flagDFS[i] = 0;
                edgeIndex[i] = -1;
            }
        }
        else // 偶数论迭代相关
        {
            for (long i = 0; i < NV; i++)
            {
                flagDFS[i] = 0;
                edgeIndex[i] = edgeStart[i + 1] - edgeStart[i];
            }
        }

        for (long i = 0; i < numUnmatchedU; i++)
        {
            long uFirst = unmatchedU[i];
            long augPathLen;
            if (iterations % 2 == 1) // odd iterations
                augPathLen = findAugPathLookahead(uFirst, flagLookahead, flagDFS, mate, G, augPath, edgeIndexLookahead, edgeIndex);
            else
                augPathLen = findAugPathLookaheadReverse(uFirst, flagLookahead, flagDFS, mate, G, augPath, edgeIndexLookahead, edgeIndex);

            // 找到增广路径增广匹配，否则记录未匹配顶点
            if (augPathLen > 0)
            {
                // augment in serial ... can be done in parallel also ...
                long u = unmatchedU[i];
                for (long k = 0; k < augPathLen; k += 2)
                {
                    mate[augPath[k]] = augPath[k + 1];
                    mate[augPath[k + 1]] = augPath[k];
                }
            }
            else
            {
                tQ[tQ_len++] = uFirst;
            }
        }

        // 如果没有未匹配顶点，或者找不到一条增广路径
        if ((tQ_len == 0) || (numUnmatchedU == tQ_len))
        {
            numUnmatchedU = 0;
            break;
        }
        // 更新未匹配顶点集合
        long *tt = tQ;
        tQ = unmatchedU;
        unmatchedU = tt;
        numUnmatchedU = tQ_len;
    }

    int card_match = 0;
    for (int i = 0; i < NV; i++)
    {
        if (mate[i] != -1)
        {
            card_match++;
        }
    }
    cout << "The cardinality of the maximum matching in current graph: " << card_match << endl;

    free(flagDFS);
    free(edgeIndex);
    free(unmatchedU);
    free(flagLookahead);
    free(edgeIndexLookahead);
    free(tQ);
    free(augPath);

    return (mate);
}