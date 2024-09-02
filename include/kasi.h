#ifndef KASI
#define KASI

#include <math.h>
#include <vector>
#include <queue>
#include <stack>
#include <ctype.h>
#include <stdio.h>
#include <chrono>
#include <ctime>
#include <unordered_set>
#include <set>
#include <random>
#include <utility>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "graphgenBP.h"
using namespace std;

#define heuristic
#define DEBUG 0
// max num nodes = 1 billion for int32
#define swap2(x, y) \
    {               \
        x = x + y;  \
        y = x - y;  \
        x = x - y;  \
    }

struct Edge
{
    int endpoint;
    int original_start;
    int original_end;

    bool operator==(const Edge &other) const
    {
        return endpoint == other.endpoint;
    }
    bool operator<(const Edge &other) const
    {
        return endpoint < other.endpoint;
    }
};

// 为了使用 std::unordered_set，我们需要特化 std::hash 结构体
namespace std
{
    template <>
    struct hash<Edge>
    {
        std::size_t operator()(const Edge &s) const
        {
            // 使用成员变量 key 来计算哈希值
            return std::hash<int>()(s.endpoint);
        }
    };
}
void printMemoryUsage();
void printheapmemory();
// 恢复图
struct Rec_Graph
{
    int n;
    std::vector<std::vector<int>>
        treeedge;
    std::vector<int> treevertex;
    std::vector<std::vector<int>> treedge_twin_1;
    std::vector<std::vector<int>> treedge_twin_2;
    std::vector<std::vector<Edge>> edge;
    int *treevertexdegree;
    int *matching;
    int *remdegree;

    // 投射与反投射
    int *project;
    int *re_project;

    // Rec_Graph() : treeedge(), treevertex(), treedge_twin_1(), treedge_twin_2(), edge()
    // {
    // }
};
typedef Rec_Graph rec_graph;

// 理论算法使用结构体
struct hash_graph
{
    std::vector<std::unordered_set<Edge>> edge; // 边表
    int *remdegree;                             // 度数表

    std::vector<std::vector<int>> treeedge;
    std::vector<int> treevertex;
    std::vector<std::vector<int>> treedge_twin_1;
    std::vector<std::vector<int>> treedge_twin_2;
    int *treevertexdegree;
    int *matching;

    // 投射与反投射
    int *project;
    int *re_project;
};

void KaSi_hash(graph *og, hash_graph *hg);

struct tree_graph
{
    std::vector<std::set<Edge>> edge; // 边表
    int *remdegree;                   // 度数表

    std::vector<std::vector<int>>
        treeedge;
    std::vector<int> treevertex;
    std::vector<std::vector<int>> treedge_twin_1;
    std::vector<std::vector<int>> treedge_twin_2;
    int *treevertexdegree;
    int *matching;

    // 投射与反投射
    int *project;
    int *re_project;
};

void KaSi_tree(graph *og, tree_graph *tg);

void KS2basic(graph *G, rec_graph *rg);
void Comp_kaSi(graph *og, rec_graph *rg);
void GMVM(graph *og, rec_graph *rg);
void BMVM(graph *og, rec_graph *rg);

// 初始化二分图
inline void init_bipartite(rec_graph *&rg, stack<int> &bucket1, stack<int> &bucket2, graph *G)
{
    rg->n = G->n;
    int rd = 0;
    Edge tmpedge;
    for (int i = 0; i < G->n; i++)
    {
        rg->remdegree[i] = G->vtx_pointer[i + 1] - G->vtx_pointer[i]; // remdegree是顶点的度数
        rg->edge[i].resize(rg->remdegree[i]);                         // edge[i]应该是表示顶点，edge[i][j]表示边，resize给顶点覆盖的边分配大小
        rd = 0;
        for (int j = G->vtx_pointer[i]; j < G->vtx_pointer[i + 1]; ++j)
        {
            tmpedge.endpoint = G->endV[j];
            tmpedge.original_start = i;
            tmpedge.original_end = G->endV[j]; // 行顶点为原始id+n
            rg->edge[i][rd++] = tmpedge;
        } // 加入边，如果度数为一或者二加入相应的桶中
        if (rg->remdegree[i] == 1)
            bucket1.push(i);
        else if (rg->remdegree[i] == 2)
            bucket2.push(i);
    }
}

inline void free_rec_graph(rec_graph *rg, int NSIZE)
{
    free(rg->matching);
    free(rg->remdegree);
    free(rg->treevertexdegree);
    rg->treevertex.clear();
    rg->treevertex.shrink_to_fit();
    for (int i = 0; i < NSIZE; ++i)
    {
        rg->treeedge[i].clear();
        rg->treeedge[i].shrink_to_fit();
        rg->treedge_twin_1[i].clear();
        rg->treedge_twin_1[i].shrink_to_fit();
        rg->treedge_twin_2[i].clear();
        rg->treedge_twin_2[i].shrink_to_fit();
        rg->edge[i].clear();
        rg->edge[i].shrink_to_fit();
    }
    rg->edge.clear();
    rg->edge.shrink_to_fit();
}

inline void free_rec_graph(hash_graph *rg, int NSIZE)
{
    free(rg->matching);
    free(rg->remdegree);
    free(rg->treevertexdegree);
    rg->treevertex.clear();
    rg->treevertex.shrink_to_fit();
    for (int i = 0; i < NSIZE; ++i)
    {
        rg->treeedge[i].clear();
        rg->treeedge[i].shrink_to_fit();
        rg->treedge_twin_1[i].clear();
        rg->treedge_twin_1[i].shrink_to_fit();
        rg->treedge_twin_2[i].clear();
        rg->treedge_twin_2[i].shrink_to_fit();
        unordered_set<Edge> emptyedgelist;
        rg->edge[i].swap(emptyedgelist);
    }
    rg->edge.clear();
    rg->edge.shrink_to_fit();
}

inline void free_rec_graph(tree_graph *rg, int NSIZE)
{
    free(rg->matching);
    free(rg->remdegree);
    free(rg->treevertexdegree);
    rg->treevertex.clear();
    rg->treevertex.shrink_to_fit();
    for (int i = 0; i < NSIZE; ++i)
    {
        rg->treeedge[i].clear();
        rg->treeedge[i].shrink_to_fit();
        rg->treedge_twin_1[i].clear();
        rg->treedge_twin_1[i].shrink_to_fit();
        rg->treedge_twin_2[i].clear();
        rg->treedge_twin_2[i].shrink_to_fit();
        set<Edge> emptyedgelist;
        rg->edge[i].swap(emptyedgelist);
    }
    rg->edge.clear();
    rg->edge.shrink_to_fit();
}
inline void remove_element(vector<Edge> &vec, int element)
{
    vec[element] = vec.back();
    vec.pop_back();
}

inline void reduce_degree(int i, int *remdegree, stack<int> &bucket1, stack<int> &bucket2)
{
    remdegree[i]--;
    if (remdegree[i] == 1)
        bucket1.push(i);
    else if (remdegree[i] == 2)
        bucket2.push(i);
}

// 移除顶点并且更新所有邻居信息
inline void basic_remove_neighborhood(int u, vector<Edge> &edge, int *remdegree, stack<int> &bucket1, stack<int> &bucket2)
{
    int examined = 0;
    int i = 0;
    while (examined < remdegree[u] - 1)
    {
        if (remdegree[edge[i].endpoint])
        {
            reduce_degree(edge[i].endpoint, remdegree, bucket1, bucket2);
            examined++;
        }
        i++;
    }
    remdegree[u] = 0;
    edge.clear();
}

// 对于pair进行hash似乎不可行（输入一个pair，返回一个hash值）
inline void reduce_the_graph(graph *og, rec_graph *rg)
{
    // 修改og,首先记录存在的顶点的数量
    rg->project = (int *)malloc(sizeof(int) * og->n);
    memset(rg->project, -1, sizeof(int) * og->n);
    rg->re_project = (int *)malloc(sizeof(int) * og->n);
    memset(rg->re_project, -1, sizeof(int) * og->n);
    int num_of_vtx = 0;
    for (int i = 0; i < og->n; i++)
    {
        if (rg->remdegree[i] != 0)
        {
            rg->project[i] = num_of_vtx;
            rg->re_project[num_of_vtx] = i;
            num_of_vtx++;
        }
    }

    // 获取行顶点的数量
    int nrows = 0;
    for (int i = og->nrows; i < og->n; i++)
    {
        if (rg->project[i] != -1)
        {
            nrows = rg->project[i];
            break;
        }
    }
    og->nrows = nrows;

    // 重新构建og
    // 处理每一个非零元素
    int actual_vtx;
    int neighbor_vtx;
    int neighbor_num;
    int examined;
    int j;
    for (int i = 0; i < num_of_vtx; i++)
    {
        // 处理每一个非零元素的邻居顶点
        actual_vtx = rg->re_project[i];
        neighbor_num = og->vtx_pointer[i];

        examined = 0, j = 0;
        while (examined < rg->remdegree[actual_vtx])
        {
            if (rg->remdegree[rg->edge[actual_vtx][j].endpoint])
            {
                og->endV[neighbor_num++] = rg->project[rg->edge[actual_vtx][j].endpoint];
                examined++;
            }
            j++;
        }
        og->vtx_pointer[i + 1] = neighbor_num;
    }

    cout << "总边数   " << neighbor_num << "顶点数   " << num_of_vtx << endl;

    og->n = num_of_vtx;
    og->m = neighbor_num;
}
// 记录匹配信息
inline void augment_recover_graph(Edge e, rec_graph *rg)
{
    int original_a, original_b;
    original_a = e.original_start;
    original_b = e.original_end;

    // 记录存在的顶点
    if (rg->treevertexdegree[original_a] == 0) // this is done to save some space
        rg->treevertex.push_back(original_a);
    rg->treevertexdegree[original_a]++;
    if (rg->treevertexdegree[original_b] == 0)
        rg->treevertex.push_back(original_b);
    rg->treevertexdegree[original_b]++;

    // 记录边        恢复图的存储仍然使用CSR类似的格式
    rg->treeedge[original_b].push_back(original_a);
    rg->treeedge[original_a].push_back(original_b);
}

inline void augment_recover_graph(Edge e, hash_graph *rg)
{
    int original_a, original_b;
    original_a = e.original_start;
    original_b = e.original_end;

    // 记录存在的顶点
    if (rg->treevertexdegree[original_a] == 0) // this is done to save some space
        rg->treevertex.push_back(original_a);
    rg->treevertexdegree[original_a]++;
    if (rg->treevertexdegree[original_b] == 0)
        rg->treevertex.push_back(original_b);
    rg->treevertexdegree[original_b]++;

    // 记录边        恢复图的存储仍然使用CSR类似的格式
    rg->treeedge[original_b].push_back(original_a);
    rg->treeedge[original_a].push_back(original_b);
}

inline void augment_recover_graph(Edge e, tree_graph *rg)
{
    int original_a, original_b;
    original_a = e.original_start;
    original_b = e.original_end;

    // 记录存在的顶点
    if (rg->treevertexdegree[original_a] == 0) // this is done to save some space
        rg->treevertex.push_back(original_a);
    rg->treevertexdegree[original_a]++;
    if (rg->treevertexdegree[original_b] == 0)
        rg->treevertex.push_back(original_b);
    rg->treevertexdegree[original_b]++;

    // 记录边        恢复图的存储仍然使用CSR类似的格式
    rg->treeedge[original_b].push_back(original_a);
    rg->treeedge[original_a].push_back(original_b);
}

inline void add_twins(Edge e1, Edge e2, rec_graph *rg)
{
    int oe11 = e1.original_start;
    int oe12 = e1.original_end;
    int oe21 = e2.original_start;
    int oe22 = e2.original_end;

    rg->treedge_twin_1[oe11].push_back(oe21);
    rg->treedge_twin_2[oe11].push_back(oe22);

    rg->treedge_twin_1[oe12].push_back(oe21);
    rg->treedge_twin_2[oe12].push_back(oe22);

    rg->treedge_twin_1[oe21].push_back(oe11);
    rg->treedge_twin_2[oe21].push_back(oe12);

    rg->treedge_twin_1[oe22].push_back(oe11);
    rg->treedge_twin_2[oe22].push_back(oe12);
}

inline void add_twins(Edge e1, Edge e2, hash_graph *rg)
{
    int oe11 = e1.original_start;
    int oe12 = e1.original_end;
    int oe21 = e2.original_start;
    int oe22 = e2.original_end;

    rg->treedge_twin_1[oe11].push_back(oe21);
    rg->treedge_twin_2[oe11].push_back(oe22);

    rg->treedge_twin_1[oe12].push_back(oe21);
    rg->treedge_twin_2[oe12].push_back(oe22);

    rg->treedge_twin_1[oe21].push_back(oe11);
    rg->treedge_twin_2[oe21].push_back(oe12);

    rg->treedge_twin_1[oe22].push_back(oe11);
    rg->treedge_twin_2[oe22].push_back(oe12);
}

inline void add_twins(Edge e1, Edge e2, tree_graph *rg)
{
    int oe11 = e1.original_start;
    int oe12 = e1.original_end;
    int oe21 = e2.original_start;
    int oe22 = e2.original_end;

    rg->treedge_twin_1[oe11].push_back(oe21);
    rg->treedge_twin_2[oe11].push_back(oe22);

    rg->treedge_twin_1[oe12].push_back(oe21);
    rg->treedge_twin_2[oe12].push_back(oe22);

    rg->treedge_twin_1[oe21].push_back(oe11);
    rg->treedge_twin_2[oe21].push_back(oe12);

    rg->treedge_twin_1[oe22].push_back(oe11);
    rg->treedge_twin_2[oe22].push_back(oe12);
}

// 恢复匹配
inline void recover_the_matching(graph *og, rec_graph *rg, long *mate)
{
    // 将缩减图的匹配转变为全局匹配
    // int temp_num1 = 0, temp_num2 = 0, temp_num3 = 0, temp_edge_num = 0;
    for (int i = 0; i < og->n; i++)
    {
        if (mate[i] != -1)
        {
            rg->matching[rg->re_project[i]] = rg->re_project[mate[i]];
            rg->matching[rg->re_project[mate[i]]] = rg->re_project[i];
            // temp_num1++;
        }
    }
    // 此处恢复的匹配与后续恢复的匹配之和不等于总匹配的原因，还有度数为一的顶点被匹配的计数没有加进来

    int incr = 0;
    // 记录度数为1的顶点和与外界相连的顶点（匹配顶点）
    stack<int> bucket1_recover;
    stack<int> bucket_matched_recover;
    for (int k = 0; k < (int)rg->treevertex.size(); ++k)
    {
        // temp_edge_num += rg->treevertexdegree[rg->treevertex[k]];
        if (rg->treevertexdegree[rg->treevertex[k]] == 1)
        {
            bucket1_recover.push(rg->treevertex[k]);
            // temp_num2++;
        }
        if (rg->matching[rg->treevertex[k]] != -1)
        {
            bucket_matched_recover.push(rg->treevertex[k]);
            // temp_num3++;
        }
    }

    // cout << temp_num1 << " " << temp_num2 << " " << temp_num3 << " " << temp_edge_num << endl;

    // int card_of_recover_matching = 0;
    // int card_of_matching1 = 0;
    // for (int i = 0; i < rg->n; i++)
    // {
    //     if (rg->matching[i] != -1)
    //     {
    //         card_of_matching1++;
    //     }
    // }
    // cout << "恢复之前的匹配为 " << card_of_matching1 << endl;

    // 将这两类顶点处理，直到为空，先处理与外界相连的顶点
    while (!bucket1_recover.empty() || !bucket_matched_recover.empty())
    {
        while (!bucket_matched_recover.empty()) // 先处理被匹配的顶点
        {
            int v = bucket_matched_recover.top();
            int mvv = rg->matching[v];
            bucket_matched_recover.pop();
            for (int i = 0; i < (int)rg->treeedge[v].size(); ++i)
            {
                int u = rg->treeedge[v][i];
                int twin_1 = rg->treedge_twin_1[v][i];
                int twin_2 = rg->treedge_twin_2[v][i];
                if (u == mvv)
                {
                    if (u > mvv)
                    {
                        rg->treevertexdegree[twin_1]--;
                        rg->treevertexdegree[twin_2]--;
                        if (rg->treevertexdegree[twin_1] == 1)
                            bucket1_recover.push(twin_1);
                        if (rg->treevertexdegree[twin_2] == 1)
                            bucket1_recover.push(twin_2);
                    }
                }
                else
                {
                    if (rg->treevertexdegree[u] > 0)
                    {
                        rg->treevertexdegree[u]--;
                        if (rg->treevertexdegree[u] == 1)
                            bucket1_recover.push(u);
                        if (rg->matching[twin_1] == -1 && rg->matching[twin_2] == -1)
                        {
                            rg->matching[twin_1] = twin_2;
                            rg->matching[twin_2] = twin_1;
                            // card_of_recover_matching++;
                            incr++;
                            bucket_matched_recover.push(twin_1);
                            bucket_matched_recover.push(twin_2);
                        }
                    }
                }
            }
            rg->treevertexdegree[v] = 0;
        }

        if (!bucket1_recover.empty()) // 然后处理度数为一的顶点
        {
            // found deg 1 vertex
            int v = bucket1_recover.top();
            bucket1_recover.pop();

            if (rg->treevertexdegree[v] > 0 && rg->matching[v] == -1)
            {
                int neigh = 0;
                while ((rg->matching[rg->treeedge[v][neigh]] != -1) || (rg->treevertexdegree[rg->treeedge[v][neigh]] == 0) || (rg->matching[rg->treedge_twin_1[v][neigh]] == rg->treedge_twin_2[v][neigh]))
                {
                    neigh++;
                    if (neigh >= (int)rg->treeedge[v].size())
                        break;
                }
                if (neigh >= (int)rg->treeedge[v].size())
                {
                    printf("this shouldnt happen!!! \n");
                    break;
                }

                // assert this
                int u = rg->treeedge[v][neigh];
                rg->matching[v] = u;
                rg->matching[u] = v;
                incr++;
                bucket_matched_recover.push(u);
                bucket_matched_recover.push(v);
                // card_of_recover_matching++;

                rg->treevertexdegree[v] = 0;
            }
        }
    }

    int card_of_matching = 0;
    for (int i = 0; i < rg->n; i++)
    {
        if (rg->matching[i] != -1)
        {
            card_of_matching++;
        }
    }
    cout << "原始图的最大基数匹配为 " << card_of_matching << endl;
    // cout << "从存储匹配中恢复匹配的有" << card_of_recover_matching << endl;

    std::ofstream of;
    of.open("TIME.txt", ios ::app);
    of << card_of_matching << "  ";
    of.close();

    free(rg->re_project);
    free(rg->project);
}

#endif // HEADER_FILE_NAME_H
