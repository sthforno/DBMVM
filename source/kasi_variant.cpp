#include <queue>
#include <stack>
#include <ctype.h>
#include <chrono>
#include <ctime>
#include <unordered_set>
#include <set>
#include <fstream>
#include <cstring>
#include <algorithm>
#include "../include/kasi.h"
#include "../include/graphgenBP.h"

using namespace std;

// max num nodes = 1 billion for int32

// 并行G里面先存储行再存储列!
// 高质量匹配启发式
void KS2basic(graph *og, rec_graph *rg)
{

    int NSIZE = og->n; // 一侧顶点的数量
    int u, w, ufu, ufw, uf1, uf2;
    stack<int> bucket1;
    stack<int> bucket2;

    Edge tmpedge;

    int remain = 0;
    int deg1count = 0;
    int deg2count = 0;
    int nnz = og->m; // 边数
    int next;
    int edges = nnz;

    rg->remdegree = (int *)malloc((NSIZE) * sizeof(int));
    rg->matching = (int *)malloc(NSIZE * sizeof(int));

    rg->edge.resize(NSIZE);

    // 恢复匹配相关
    rg->treeedge.resize(NSIZE);
    rg->treedge_twin_1.resize(NSIZE);
    rg->treedge_twin_2.resize(NSIZE);
    rg->treevertexdegree = (int *)calloc(NSIZE, sizeof(int));

    // 初始化匹配
    memset(rg->matching, -1, sizeof(int) * NSIZE);
    int *mergearray;
    mergearray = (int *)calloc(NSIZE, sizeof(int));
    // printf("use merge array \n");

    init_bipartite(rg, bucket1, bucket2, og);

    int merge_operations = 0;
    // 主循环，直到处理完所有的边? 并且度数一和度数二的桶都为空
    while (!bucket1.empty() || !bucket2.empty()) //
    {
        // 度数一的顶点
        while (!bucket1.empty())
        {
            int v = bucket1.top();
            bucket1.pop();
            if (rg->remdegree[v] == 1)
            {
                deg1count++;
                int point = 0;
                // 找到存在的邻居顶点
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                u = rg->edge[v][point].endpoint;
                // 将该边的原始边加入匹配
                rg->matching[rg->edge[v][point].original_start] = rg->edge[v][point].original_end;
                rg->matching[rg->edge[v][point].original_end] = rg->edge[v][point].original_start;

                // 是被合并边，原始端点不一定被移除了。
                rg->remdegree[v] = 0;
                basic_remove_neighborhood(u, rg->edge[u], rg->remdegree, bucket1, bucket2);
                rg->edge[v].clear();
            }
        }

        if (!bucket2.empty())
        {
            int v = bucket2.top();
            bucket2.pop();
            if (rg->remdegree[v] == 2)
            {
                merge_operations++;
                // 获取两个端点及相关的边
                int point = 0;
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                int u = rg->edge[v][point].endpoint;
                Edge e1 = rg->edge[v][point];
                point++;
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                int w = rg->edge[v][point].endpoint;
                Edge e2 = rg->edge[v][point];

                // 将度数较大的设置为u
                if (rg->remdegree[w] > rg->remdegree[u])
                    swap2(u, w);

                // // 构建recover图，用于恢复匹配
                augment_recover_graph(e1, rg);
                augment_recover_graph(e2, rg);
                add_twins(e1, e2, rg);

                rg->remdegree[v] = 0; // 删除顶点v

                // 处理顶点u，将存在的顶点进行记录，移除被删除的边
                int i = 0;
                int upper_limit = rg->edge[u].size() - 1;
                for (int i = upper_limit; i >= 0; --i)
                {
                    if (rg->remdegree[rg->edge[u][i].endpoint])
                    {
                        mergearray[rg->edge[u][i].endpoint] = u;
                    }
                    else
                        remove_element(rg->edge[u], i);
                }

                // 处理顶点w，对于w的邻居，如果是非并行边，将其进行两次插入（两个顶点）
                // 如果是并行边，更新度数，并且判断是否需要重新处理
                rg->remdegree[w]--;
                int examined = 0;
                i = 0;
                while (examined < rg->remdegree[w])
                {
                    if (rg->remdegree[rg->edge[w][i].endpoint])
                    {
                        if (mergearray[rg->edge[w][i].endpoint] != u)
                        {
                            rg->edge[u].push_back(rg->edge[w][i]);
                            tmpedge.endpoint = u;
                            tmpedge.original_start = rg->edge[w][i].original_end;
                            tmpedge.original_end = rg->edge[w][i].original_start;
                            rg->edge[rg->edge[w][i].endpoint].push_back(tmpedge);
                            rg->remdegree[u]++;
                        }
                        else
                        {
                            reduce_degree(rg->edge[w][i].endpoint, rg->remdegree, bucket1, bucket2);
                        }
                        examined++;
                    }
                    i++;
                }

                // 最终结果，删除了w和v顶点
                rg->remdegree[w] = 0;
                deg2count++;
                reduce_degree(u, rg->remdegree, bucket1, bucket2);
                rg->edge[v].clear();
                rg->edge[w].clear();
            }
        }
    }
    auto end_kasi = std::chrono::system_clock::now();

    int num_of_vtx = 0;
    for (int i = 0; i < og->n; i++)
    {
        if (rg->remdegree[i] > 0)
        {
            num_of_vtx++;
        }
    }
    cout << "剩下的顶点数量为" << num_of_vtx << endl;
    cout << "获取的匹配基数为" << deg1count + deg2count << endl;

    // 记录顶点数，剩余顶点数，合并次数
    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // // of << og->n << "    " << num_of_vtx << "    " << deg2count << "    ";
    // of << merge_operations << "  ";
    // of.close();
    // printMemoryUsage();
    // printheapmemory();
    free(mergearray);
}

#define CACHE_SIZE 5
struct cache_vertices
{
    std::vector<std::unordered_set<int>> high_degree_vertices_cache;
    std::vector<int> high_degree_vertices_cache_id;
    bool is_high_degree_vertices_cache_empty = true;
};

void Cache_high_degree_vertices(graph *og, rec_graph *rg, cache_vertices *ca, int v, int *frequency)
{
    if (ca->is_high_degree_vertices_cache_empty = true)
    {
        for (int i = 0; i < CACHE_SIZE; i++)
        {
            if (ca->high_degree_vertices_cache_id[i] == -1)
            {
                ca->high_degree_vertices_cache_id[i] = v;
                ca->high_degree_vertices_cache[i].clear();
                for (int j = 0; j < rg->edge[v].size(); j++)
                {
                    if (rg->remdegree[rg->edge[v][j].endpoint] > 0)
                        ca->high_degree_vertices_cache[i].insert(rg->edge[v][j].endpoint);
                }
                break;
            }
            if (i == CACHE_SIZE - 1)
            {
                ca->is_high_degree_vertices_cache_empty = false;
            }
        }
    }
    else // 从cache中移除数据需要更新边表
    {
        for (int i = 0; i < CACHE_SIZE; i++)
        {
            if (frequency[ca->high_degree_vertices_cache_id[i]] < frequency[v])
            {
                ca->high_degree_vertices_cache_id[i] = v;
                ca->high_degree_vertices_cache[i].clear();
                for (int j = 0; j < rg->edge[v].size(); j++)
                {
                    if (rg->remdegree[rg->edge[v][j].endpoint] > 0)
                        ca->high_degree_vertices_cache[i].insert(rg->edge[v][j].endpoint);
                }
                break;
            }
        }
    }
}

int in_cache(int v, cache_vertices *ca)
{
    int i;
    for (i = 0; i < CACHE_SIZE; i++)
    {
        if (ca->high_degree_vertices_cache_id[i] == v)
        {
            return i;
        }
    }
    return CACHE_SIZE;
}

// 只有在合并顶点时会反复访问高度顶点所有邻居
void Cache_KaSi(graph *og, rec_graph *rg)
{

    int NSIZE = og->n;
    int u, w, ufu, ufw, uf1, uf2;
    stack<int> bucket1;
    stack<int> bucket2;

    Edge tmpedge;

    int remain = 0;
    int deg1count = 0;
    int deg2count = 0;
    int nnz = og->m; // 边数
    int next;
    int edges = nnz;

    int degree_threshold = log(NSIZE);

    rg->remdegree = (int *)malloc((NSIZE) * sizeof(int));
    rg->matching = (int *)malloc(NSIZE * sizeof(int));

    rg->edge.resize(NSIZE);

    // 恢复匹配相关
    rg->treeedge.resize(NSIZE);
    rg->treedge_twin_1.resize(NSIZE);
    rg->treedge_twin_2.resize(NSIZE);
    rg->treevertexdegree = (int *)calloc(NSIZE, sizeof(int));

    // 初始化匹配
    memset(rg->matching, -1, sizeof(int) * NSIZE);
    int *mergearray;
    mergearray = (int *)calloc(NSIZE, sizeof(int));
    // printf("use merge array \n");

    init_bipartite(rg, bucket1, bucket2, og);

    // cache频繁处理的高度顶点
    cache_vertices ca;
    ca.high_degree_vertices_cache.resize(CACHE_SIZE);
    ca.high_degree_vertices_cache_id.resize(CACHE_SIZE);
    ca.high_degree_vertices_cache_id.assign(CACHE_SIZE, -1);

    int *vertices_frequency = (int *)malloc(sizeof(int) * NSIZE);
    memset(vertices_frequency, 0, sizeof(int) * NSIZE);

    int merge_operations = 0;
    // 主循环，直到处理完所有的边? 并且度数一和度数二的桶都为空
    while (!bucket1.empty() || !bucket2.empty()) //
    {
        // 度数一的顶点
        while (!bucket1.empty())
        {
            int v = bucket1.top();
            bucket1.pop();
            if (rg->remdegree[v] == 1)
            {
                deg1count++;
                int point = 0;
                // 找到存在的邻居顶点
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                u = rg->edge[v][point].endpoint;
                // 将该边的原始边加入匹配
                rg->matching[rg->edge[v][point].original_start] = rg->edge[v][point].original_end;
                rg->matching[rg->edge[v][point].original_end] = rg->edge[v][point].original_start;

                // 是被合并边，原始端点不一定被移除了。
                rg->remdegree[v] = 0;
                basic_remove_neighborhood(u, rg->edge[u], rg->remdegree, bucket1, bucket2);
                rg->edge[v].clear();
            }
        }

        if (!bucket2.empty())
        {
            int v = bucket2.top();
            bucket2.pop();
            if (rg->remdegree[v] == 2)
            {
                merge_operations++;
                // 获取两个端点及相关的边
                int point = 0;
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                int u = rg->edge[v][point].endpoint;
                Edge e1 = rg->edge[v][point];
                point++;
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                int w = rg->edge[v][point].endpoint;
                Edge e2 = rg->edge[v][point];

                // 将度数较大的设置为u
                if (rg->remdegree[w] > rg->remdegree[u])
                    swap2(u, w);

                vertices_frequency[u]++;
                vertices_frequency[w]++;
                // 如果cache没满，直接找空插入
                if (rg->remdegree[u] > degree_threshold)
                {
                    if (in_cache(u, &ca) == CACHE_SIZE)
                        Cache_high_degree_vertices(og, rg, &ca, u, vertices_frequency);
                }
                if (rg->remdegree[w] > degree_threshold)
                {
                    if (in_cache(w, &ca) == CACHE_SIZE)
                        Cache_high_degree_vertices(og, rg, &ca, w, vertices_frequency);
                }

                // // 构建recover图，用于恢复匹配
                augment_recover_graph(e1, rg);
                augment_recover_graph(e2, rg);
                add_twins(e1, e2, rg);

                rg->remdegree[v] = 0; // 删除顶点v

                // 处理顶点u，将存在的顶点进行记录，移除被删除的边
                int position = in_cache(u, &ca);
                if (position < CACHE_SIZE)
                {
                    // 对于u中已经被移除的元素怎么处理？没法处理，只能够在每次遍历的时候加上判断语句是否度数为0
                    // 处理顶点w
                    rg->remdegree[w]--;
                    int examined = 0;
                    int i = 0;
                    while (examined < rg->remdegree[w])
                    {
                        if (rg->remdegree[rg->edge[w][i].endpoint])
                        {
                            // 如果在cache中能够找到元素
                            if (ca.high_degree_vertices_cache[position].find(rg->edge[w][i].endpoint) != ca.high_degree_vertices_cache[position].end())
                            {
                                reduce_degree(rg->edge[w][i].endpoint, rg->remdegree, bucket1, bucket2);
                            }
                            else // 找不到元素，同时插入cache和原始的邻接表
                            {
                                ca.high_degree_vertices_cache[position].insert(rg->edge[w][i].endpoint);
                                rg->edge[u].push_back(rg->edge[w][i]);
                                tmpedge.endpoint = u;
                                tmpedge.original_start = rg->edge[w][i].original_end;
                                tmpedge.original_end = rg->edge[w][i].original_start;
                                rg->edge[rg->edge[w][i].endpoint].push_back(tmpedge);
                                rg->remdegree[u]++;
                            }
                            examined++;
                        }
                        i++;
                    }

                    // 最终结果，删除了w和v顶点
                    rg->remdegree[w] = 0;
                    deg2count++;
                    reduce_degree(u, rg->remdegree, bucket1, bucket2);
                    rg->edge[v].clear();
                    rg->edge[w].clear();
                }
                else // 如果不在cache中，正常更新
                {
                    int upper_limit = rg->edge[u].size() - 1;
                    for (int i = upper_limit; i >= 0; --i)
                    {
                        if (rg->remdegree[rg->edge[u][i].endpoint])
                        {
                            mergearray[rg->edge[u][i].endpoint] = u;
                        }
                        else
                            remove_element(rg->edge[u], i);
                    }

                    // 处理顶点w，对于w的邻居，如果是非并行边，将其进行两次插入（两个顶点）
                    // 如果是并行边，更新度数，并且判断是否需要重新处理
                    rg->remdegree[w]--;
                    int examined = 0;
                    int i = 0;
                    while (examined < rg->remdegree[w])
                    {
                        if (rg->remdegree[rg->edge[w][i].endpoint])
                        {
                            if (mergearray[rg->edge[w][i].endpoint] != u)
                            {
                                rg->edge[u].push_back(rg->edge[w][i]);
                                tmpedge.endpoint = u;
                                tmpedge.original_start = rg->edge[w][i].original_end;
                                tmpedge.original_end = rg->edge[w][i].original_start;
                                rg->edge[rg->edge[w][i].endpoint].push_back(tmpedge);
                                rg->remdegree[u]++;
                            }
                            else
                            {
                                reduce_degree(rg->edge[w][i].endpoint, rg->remdegree, bucket1, bucket2);
                            }
                            examined++;
                        }
                        i++;
                    }

                    // 最终结果，删除了w和v顶点
                    rg->remdegree[w] = 0;
                    deg2count++;
                    reduce_degree(u, rg->remdegree, bucket1, bucket2);
                    rg->edge[v].clear();
                    rg->edge[w].clear();
                }
            }
        }
    }
    auto end_kasi = std::chrono::system_clock::now();

    int num_of_vtx = 0;
    for (int i = 0; i < og->n; i++)
    {
        if (rg->remdegree[i] > 0)
        {
            num_of_vtx++;
        }
    }
    cout << "剩下的顶点数量为" << num_of_vtx << endl;
    cout << "获取的匹配基数为" << deg1count + deg2count << endl;

    // 记录顶点数，剩余顶点数，合并次数
    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // // of << og->n << "    " << num_of_vtx << "    " << deg2count << "    ";
    // of << merge_operations << "  ";
    // of.close();
    // printMemoryUsage();
    printheapmemory();
    free(mergearray);
    free(vertices_frequency);
}

//==================================================================================================================================================//
//==========================================理论算法一：使用hash表存储图结构===========================================================================//
//==================================================================================================================================================//

// 初始化hash图
void init_hashgraph(graph *og, hash_graph *hg, stack<int> &bucket1, stack<int> &bucket2)
{
    int rd = 0;
    Edge tmpedge;
    hg->remdegree = (int *)malloc((og->n) * sizeof(int));
    hg->matching = (int *)malloc(og->n * sizeof(int));
    memset(hg->matching, -1, sizeof(int) * og->n);
    hg->edge.resize(og->n);
    for (int i = 0; i < og->n; i++)
    {
        hg->remdegree[i] = og->vtx_pointer[i + 1] - og->vtx_pointer[i]; // 先获取度数，再存入边
        rd = 0;
        for (int j = og->vtx_pointer[i]; j < og->vtx_pointer[i + 1]; ++j)
        {
            tmpedge.endpoint = og->endV[j];
            tmpedge.original_start = i;
            tmpedge.original_end = og->endV[j]; // 行顶点为原始id+n
            hg->edge[i].insert(tmpedge);        // 存入边
        }
        if (hg->remdegree[i] == 1)
            bucket1.push(i);
        else if (hg->remdegree[i] == 2)
            bucket2.push(i);
    }

    // 恢复匹配相关
    hg->treeedge.resize(og->n);
    hg->treedge_twin_1.resize(og->n);
    hg->treedge_twin_2.resize(og->n);
    hg->treevertexdegree = (int *)calloc(og->n, sizeof(int));
}

void KaSi_hash(graph *og, hash_graph *hg)
{
    // 两个存储可以处理顶点的桶
    stack<int> bucket1;
    stack<int> bucket2;
    // 两类顶点操作计数
    int deg1count = 0;
    int deg2count = 0;

    // 对于hash图进行初始化

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // auto start_kasi = std::chrono::system_clock::now();
    init_hashgraph(og, hg, bucket1, bucket2);
    // auto end_kasi = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // of << elapsed_seconds_kasi.count() << "    ";
    // of.close();
    // 主循环，直到处理完所有的边? 并且度数一和度数二的桶都为空
    while (!bucket1.empty() || !bucket2.empty()) //
    {
        // 度数一的顶点
        while (!bucket1.empty())
        {
            int v = bucket1.top();
            bucket1.pop();
            if (hg->remdegree[v] == 1)
            {
                deg1count++;

                int u = -1;
                Edge record_edge;
                for (const auto &edge : hg->edge[v])
                {
                    if (hg->remdegree[edge.endpoint])
                    {
                        u = edge.endpoint;
                        record_edge = edge;
                        break;
                    }
                }

                hg->matching[record_edge.original_start] = record_edge.original_end;
                hg->matching[record_edge.original_end] = record_edge.original_start;

                hg->remdegree[v] = 0; // 同时移除这两个顶点

                // 更新所有邻居顶点
                for (const auto &neighbor_edge : hg->edge[u])
                {
                    if (hg->remdegree[neighbor_edge.endpoint])
                    {
                        // 处理u的邻居顶点，移除u
                        auto it = hg->edge[neighbor_edge.endpoint].find(record_edge);
                        if (it != hg->edge[neighbor_edge.endpoint].end())
                        {
                            hg->edge[neighbor_edge.endpoint].erase(it); // 找到了元素并移除
                        }
                        hg->remdegree[neighbor_edge.endpoint]--;
                        if (hg->remdegree[neighbor_edge.endpoint] == 1)
                            bucket1.push(neighbor_edge.endpoint);
                        else if (hg->remdegree[neighbor_edge.endpoint] == 2)
                            bucket2.push(neighbor_edge.endpoint);
                    }
                }

                hg->remdegree[u] = 0;
                hg->edge[u].clear();
                hg->edge[v].clear();
            }
        }

        if (!bucket2.empty())
        {
            int v = bucket2.top();
            bucket2.pop();
            if (hg->remdegree[v] == 2)
            {
                int temp_vtx_que[2];
                int temp_que_end = 0;
                vector<Edge> temp_edge_que;
                temp_edge_que.clear();
                for (const auto &edge : hg->edge[v])
                {
                    if (hg->remdegree[edge.endpoint])
                    {
                        temp_vtx_que[temp_que_end++] = edge.endpoint;
                        temp_edge_que.push_back(edge);
                        if (temp_que_end == 2)
                            break;
                    }
                }

                hg->remdegree[v] = 0;
                hg->edge[v].clear();
                augment_recover_graph(temp_edge_que[0], hg);
                augment_recover_graph(temp_edge_que[1], hg);
                add_twins(temp_edge_que[0], temp_edge_que[1], hg);

                int u = temp_vtx_que[0];
                int w = temp_vtx_que[1];

                // 不需要记录顶点u的连接状态，直接处理w中的每一个顶点即可

                // 处理顶点w，对于w的邻居，如果是非并行边，将其进行两次插入（两个顶点）
                // 如果是并行边，更新度数，并且判断是否需要重新处理
                for (const auto &neighbor_edge : hg->edge[w])
                {
                    if (hg->remdegree[neighbor_edge.endpoint])
                    {
                        auto item = hg->edge[u].find(neighbor_edge);
                        if (item != hg->edge[u].end()) // 存在于u邻居集合中，更新邻居顶点的度数
                        {
                            hg->remdegree[neighbor_edge.endpoint]--;
                            hg->edge[neighbor_edge.endpoint].erase(temp_edge_que[1]);

                            if (hg->remdegree[neighbor_edge.endpoint] == 1)
                                bucket1.push(neighbor_edge.endpoint);
                            else if (hg->remdegree[neighbor_edge.endpoint] == 2)
                                bucket2.push(neighbor_edge.endpoint);
                        }
                        else // 不存在与u的邻居集合中，更新u邻居列表
                        {
                            // 连接cb
                            hg->edge[u].insert(neighbor_edge);

                            // 连接bc
                            Edge tmpedge;
                            tmpedge.endpoint = u;
                            tmpedge.original_start = neighbor_edge.original_end;
                            tmpedge.original_end = neighbor_edge.original_start;
                            hg->edge[neighbor_edge.endpoint].insert(tmpedge);
                            // 移除ba
                            hg->edge[neighbor_edge.endpoint].erase(temp_edge_que[1]); // 断开原来的边并且重新连接
                            hg->remdegree[u]++;
                        }
                    }
                }

                // 最终结果，删除了w和v顶点
                hg->remdegree[w] = 0;
                hg->edge[w].clear();
                deg2count++;
                hg->remdegree[u]--;
                if (hg->remdegree[u] == 1)
                    bucket1.push(u);
                else if (hg->remdegree[u] == 2)
                    bucket2.push(u);
            }
        }
    }

    // printMemoryUsage();
    // printheapmemory();

    // 精确算法
    // 对于度数为1的顶点，直接设置相关的边的顶点的度数为0表示删除，并且更新相关顶点的度数没有实际更新边表。
    // 对于度数为2的顶点，修改其中一个顶点的连接状态，删除另外两个顶点。        通过head记录被合并到的顶点
    // 度数为0的顶点表示被删除

    // cout << "获取的匹配基数为" << deg1count + deg2count << endl;
}

//==================================================================================================================================================//
//==========================================理论算法二：使用二叉树存储图结构===========================================================================//
//==================================================================================================================================================//

// // 初始化hash图
void init_treegraph(graph *og, tree_graph *tg, stack<int> &bucket1, stack<int> &bucket2)
{
    int rd = 0;
    Edge tmpedge;
    tg->remdegree = (int *)malloc((og->n) * sizeof(int));
    tg->matching = (int *)malloc(og->n * sizeof(int));
    memset(tg->matching, -1, sizeof(int) * og->n);
    tg->edge.resize(og->n);

    for (int i = 0; i < og->n; i++)
    {
        tg->remdegree[i] = og->vtx_pointer[i + 1] - og->vtx_pointer[i]; // 先获取度数，再存入边
        rd = 0;
        for (int j = og->vtx_pointer[i]; j < og->vtx_pointer[i + 1]; ++j)
        {
            tmpedge.endpoint = og->endV[j];
            tmpedge.original_start = i;
            tmpedge.original_end = og->endV[j]; // 行顶点为原始id+n
            tg->edge[i].insert(tmpedge);        // 存入边
        }
        if (tg->remdegree[i] == 1)
            bucket1.push(i);
        else if (tg->remdegree[i] == 2)
            bucket2.push(i);
    }
    tg->treeedge.resize(og->n);
    tg->treedge_twin_1.resize(og->n);
    tg->treedge_twin_2.resize(og->n);
    tg->treevertexdegree = (int *)calloc(og->n, sizeof(int));
}

void KaSi_tree(graph *og, tree_graph *tg)
{
    // 两个存储可以处理顶点的桶
    stack<int> bucket1;
    stack<int> bucket2;
    // 两类顶点操作计数
    int deg1count = 0;
    int deg2count = 0;

    // 对于hash图进行初始化

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // auto start_kasi = std::chrono::system_clock::now();
    init_treegraph(og, tg, bucket1, bucket2);
    // auto end_kasi = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // of << elapsed_seconds_kasi.count() << endl;
    // of.close();

    // 主循环，直到处理完所有的边? 并且度数一和度数二的桶都为空
    while (!bucket1.empty() || !bucket2.empty()) //
    {
        // 度数一的顶点
        while (!bucket1.empty())
        {
            int v = bucket1.top();
            bucket1.pop();
            if (tg->remdegree[v] == 1)
            {
                deg1count++;

                int u = -1;
                Edge record_edge;
                for (const auto &edge : tg->edge[v])
                {
                    if (tg->remdegree[edge.endpoint])
                    {
                        u = edge.endpoint;
                        record_edge = edge;
                        break;
                    }
                }
                tg->matching[record_edge.original_start] = record_edge.original_end;
                tg->matching[record_edge.original_end] = record_edge.original_start;

                tg->remdegree[v] = 0; // 同时移除这两个顶点

                for (const auto &neighbor_edge : tg->edge[u])
                {
                    if (tg->remdegree[neighbor_edge.endpoint])
                    {
                        auto it = tg->edge[neighbor_edge.endpoint].find(record_edge);
                        if (it != tg->edge[neighbor_edge.endpoint].end())
                        {
                            tg->edge[neighbor_edge.endpoint].erase(it); // 找到了元素并移除
                        }
                        tg->remdegree[neighbor_edge.endpoint]--;
                        if (tg->remdegree[neighbor_edge.endpoint] == 1)
                            bucket1.push(neighbor_edge.endpoint);
                        else if (tg->remdegree[neighbor_edge.endpoint] == 2)
                            bucket2.push(neighbor_edge.endpoint);
                    }
                }
                tg->remdegree[u] = 0;
                tg->edge[u].clear();
                tg->edge[v].clear();
            }
        }

        if (!bucket2.empty())
        {
            int v = bucket2.top();
            bucket2.pop();
            if (tg->remdegree[v] == 2)
            {
                int temp_vtx_que[2];
                int temp_que_end = 0;
                vector<Edge> temp_edge_que;
                temp_edge_que.clear();
                for (const auto &edge : tg->edge[v])
                {
                    if (tg->remdegree[edge.endpoint])
                    {
                        temp_vtx_que[temp_que_end++] = edge.endpoint;
                        temp_edge_que.push_back(edge);
                        if (temp_que_end == 2)
                            break;
                    }
                }

                tg->remdegree[v] = 0;
                tg->edge[v].clear();
                augment_recover_graph(temp_edge_que[0], tg);
                augment_recover_graph(temp_edge_que[1], tg);
                add_twins(temp_edge_que[0], temp_edge_que[1], tg);

                int u = temp_vtx_que[0];
                int w = temp_vtx_que[1];

                // 不需要记录顶点u的连接状态，直接处理w中的每一个顶点即可

                // 处理顶点w，对于w的邻居，如果是非并行边，将其进行两次插入（两个顶点）
                // 如果是并行边，更新度数，并且判断是否需要重新处理
                for (const auto &neighbor_edge : tg->edge[w])
                {
                    if (tg->remdegree[neighbor_edge.endpoint])
                    {
                        auto item = tg->edge[u].find(neighbor_edge);
                        if (item != tg->edge[u].end()) // 存在于u邻居集合中，更新邻居顶点的度数
                        {
                            tg->remdegree[neighbor_edge.endpoint]--;
                            tg->edge[neighbor_edge.endpoint].erase(temp_edge_que[1]);

                            if (tg->remdegree[neighbor_edge.endpoint] == 1)
                                bucket1.push(neighbor_edge.endpoint);
                            else if (tg->remdegree[neighbor_edge.endpoint] == 2)
                                bucket2.push(neighbor_edge.endpoint);
                        }
                        else // 不存在与u的邻居集合中，更新u邻居列表
                        {
                            // 连接cb
                            tg->edge[u].insert(neighbor_edge);

                            // 连接bc
                            Edge tmpedge;
                            tmpedge.endpoint = u;
                            tmpedge.original_start = neighbor_edge.original_end;
                            tmpedge.original_end = neighbor_edge.original_start;
                            tg->edge[neighbor_edge.endpoint].insert(tmpedge);
                            // 移除ba
                            tg->edge[neighbor_edge.endpoint].erase(temp_edge_que[1]); // 断开原来的边并且重新连接
                            tg->remdegree[u]++;
                        }
                    }
                }

                // 最终结果，删除了w和v顶点
                tg->remdegree[w] = 0;
                tg->edge[w].clear();
                deg2count++;
                tg->remdegree[u]--;
                if (tg->remdegree[u] == 1)
                    bucket1.push(u);
                else if (tg->remdegree[u] == 2)
                    bucket2.push(u);
            }
        }
    }

    // printMemoryUsage();
    // printheapmemory();

    // 精确算法
    // 对于度数为1的顶点，直接设置相关的边的顶点的度数为0表示删除，并且更新相关顶点的度数没有实际更新边表。
    // 对于度数为2的顶点，修改其中一个顶点的连接状态，删除另外两个顶点。        通过head记录被合并到的顶点
    // 度数为0的顶点表示被删除

    // cout << "获取的匹配基数为" << deg1count + deg2count << endl;
}
