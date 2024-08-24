#include <queue>
#include <stack>
#include <ctype.h>
#include <chrono>
#include <ctime>
#include <unordered_set>
#include <set>
#include <fstream>
#include <cstring>
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

    // 精确算法
    // 对于度数为1的顶点，直接设置相关的边的顶点的度数为0表示删除，并且更新相关顶点的度数没有实际更新边表。
    // 对于度数为2的顶点，修改其中一个顶点的连接状态，删除另外两个顶点。        通过head记录被合并到的顶点
    // 度数为0的顶点表示被删除
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
    free(mergearray);
}

void Comp_kaSi(graph *og, rec_graph *rg)
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

    int *mergearray = (int *)calloc(NSIZE, sizeof(int));
    // printf("use merge array \n");

    init_bipartite(rg, bucket1, bucket2, og);

    int merge_operations = 0;

    int *comp_record = (int *)malloc(sizeof(int) * NSIZE);
    memset(comp_record, -1, sizeof(int) * NSIZE);

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
            comp_record[v] = v;
            bucket2.pop();
            if (rg->remdegree[v] == 2)
            {
                merge_operations++;
                // 获取两个端点及相关的边
                int point = 0;
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                int u = rg->edge[v][point].endpoint;
                comp_record[u] = v;
                Edge e1 = rg->edge[v][point];
                point++;
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                int w = rg->edge[v][point].endpoint;
                comp_record[w] = v;
                Edge e2 = rg->edge[v][point];

                augment_recover_graph(e1, rg);
                augment_recover_graph(e2, rg);
                add_twins(e1, e2, rg);

                // 寻找边界顶点
                int temp_vtx, temp_vtx1;
                while (rg->remdegree[u] == 2)
                {
                    // 获取u的另外一个邻居
                    int ptr = 0;
                    // 如果被处理过了，继续处理，如果度数为0，继续处理
                    while ((!rg->remdegree[rg->edge[u][ptr].endpoint] || comp_record[rg->edge[u][ptr].endpoint] == v))
                    {
                        ptr++;
                        if (ptr == rg->edge[u].size())
                        {
                            break;
                        }
                    }

                    if (ptr == rg->edge[u].size())
                    {
                        break;
                    }
                    temp_vtx = rg->edge[u][ptr].endpoint;
                    Edge e1 = rg->edge[u][ptr];
                    comp_record[temp_vtx] = v;
                    if (rg->remdegree[temp_vtx] != 2)
                    {
                        break;
                    }

                    // 判断u的另外一个邻居x的度数是否为2，为2继续进行查找，否则，直接以该顶点作为端点
                    ptr = 0;
                    while ((!rg->remdegree[rg->edge[temp_vtx][ptr].endpoint] || comp_record[rg->edge[temp_vtx][ptr].endpoint] == v))
                    {
                        ptr++;
                        if (ptr == rg->edge[temp_vtx].size())
                        {
                            break;
                        }
                    }

                    if (ptr == rg->edge[temp_vtx].size())
                    {
                        break;
                    }

                    temp_vtx1 = rg->edge[temp_vtx][ptr].endpoint;
                    Edge e2 = rg->edge[temp_vtx][ptr];
                    comp_record[temp_vtx1] = v;

                    // 移除路径中间的两个顶点
                    rg->remdegree[u] = 0;
                    rg->remdegree[temp_vtx] = 0;

                    augment_recover_graph(e1, rg);
                    augment_recover_graph(e2, rg);
                    add_twins(e1, e2, rg);

                    // 更新顶点u
                    u = temp_vtx1;
                }

                while (rg->remdegree[w] == 2)
                {
                    int ptr = 0;
                    // 如果被处理过了，继续处理，如果度数为0，继续处理
                    while ((!rg->remdegree[rg->edge[w][ptr].endpoint] || comp_record[rg->edge[w][ptr].endpoint] == v))
                    {
                        ptr++;
                        if (ptr == rg->edge[w].size())
                        {
                            break;
                        }
                    }

                    if (ptr == rg->edge[w].size())
                    {
                        break;
                    }

                    temp_vtx = rg->edge[w][ptr].endpoint;
                    Edge e1 = rg->edge[w][ptr];
                    comp_record[temp_vtx] = v;
                    if (rg->remdegree[temp_vtx] != 2)
                    {
                        break;
                    }

                    // 判断u的另外一个邻居x的度数是否为2，为2继续进行查找，否则，直接以该顶点作为端点
                    ptr = 0;
                    while ((!rg->remdegree[rg->edge[temp_vtx][ptr].endpoint] || comp_record[rg->edge[temp_vtx][ptr].endpoint] == v))
                    {
                        ptr++;
                        if (ptr == rg->edge[temp_vtx].size())
                        {
                            break;
                        }
                    }

                    if (ptr == rg->edge[temp_vtx].size())
                    {
                        break;
                    }

                    temp_vtx1 = rg->edge[temp_vtx][ptr].endpoint;
                    Edge e2 = rg->edge[temp_vtx][ptr];
                    comp_record[temp_vtx1] = v;

                    // 移除路径中间的两个顶点
                    rg->remdegree[w] = 0;
                    rg->remdegree[temp_vtx] = 0;

                    augment_recover_graph(e1, rg);
                    augment_recover_graph(e2, rg);
                    add_twins(e1, e2, rg);

                    // 更新顶点u
                    w = temp_vtx1;
                }

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

    // 精确算法
    // 对于度数为1的顶点，直接设置相关的边的顶点的度数为0表示删除，并且更新相关顶点的度数没有实际更新边表。
    // 对于度数为2的顶点，修改其中一个顶点的连接状态，删除另外两个顶点。        通过head记录被合并到的顶点
    // 度数为0的顶点表示被删除
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
    free(mergearray);
    free(comp_record);
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

    // 精确算法
    // 对于度数为1的顶点，直接设置相关的边的顶点的度数为0表示删除，并且更新相关顶点的度数没有实际更新边表。
    // 对于度数为2的顶点，修改其中一个顶点的连接状态，删除另外两个顶点。        通过head记录被合并到的顶点
    // 度数为0的顶点表示被删除

    cout << "获取的匹配基数为" << deg1count + deg2count << endl;
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

    // 精确算法
    // 对于度数为1的顶点，直接设置相关的边的顶点的度数为0表示删除，并且更新相关顶点的度数没有实际更新边表。
    // 对于度数为2的顶点，修改其中一个顶点的连接状态，删除另外两个顶点。        通过head记录被合并到的顶点
    // 度数为0的顶点表示被删除

    cout << "获取的匹配基数为" << deg1count + deg2count << endl;
}

void BMVM(graph *og, rec_graph *rg)
{

    int NSIZE = og->n; // 一侧顶点的数量
    int u, w;
    stack<int> bucket1;
    stack<int> bucket2;
    stack<int> bucket3;

    Edge tmpedge;
    int nnz = og->m;

    rg->remdegree = (int *)malloc(NSIZE * sizeof(int));
    rg->edge.resize(NSIZE);
    rg->matching = (int *)malloc(NSIZE * sizeof(int));

    rg->treeedge.resize(NSIZE);
    rg->treedge_twin_1.resize(NSIZE);
    rg->treedge_twin_2.resize(NSIZE);
    rg->treevertexdegree = (int *)calloc(NSIZE, sizeof(int));

    // 初始化匹配
    memset(rg->matching, -1, sizeof(int) * NSIZE);
    int *mergearray;
    mergearray = (int *)calloc(NSIZE, sizeof(int));

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // auto start_kasi = std::chrono::system_clock::now();
    init_bipartite(rg, bucket1, bucket2, og);
    // auto end_kasi = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // of << elapsed_seconds_kasi.count() << "    ";
    // of.close();

    // 存储搜索图，对于搜索图中顶点使用bfs，因此使用记录
    int *row_que = (int *)malloc(sizeof(int) * og->nrows);
    int *col_que = (int *)malloc(sizeof(int) * og->nrows);
    int row_que_start, row_que_end, col_que_start, col_que_end;
    // 记录访问顶点               1：记录是否被处理过（记录起始顶点），2：记录被访问次数，3：记录是否已经被加入子图
    int **visited = (int **)malloc(sizeof(int *) * 3);
    visited[0] = (int *)malloc(sizeof(int) * NSIZE);
    visited[1] = (int *)malloc(sizeof(int) * NSIZE);
    visited[2] = (int *)malloc(sizeof(int) * NSIZE);
    memset(visited[0], -1, sizeof(int) * NSIZE);
    memset(visited[1], -1, sizeof(int) * NSIZE);
    memset(visited[2], -1, sizeof(int) * NSIZE);

    int temp_num = 0;
    int deg1count = 0;
    int deg2count = 0;

    // 按照轮次处理(两种轮次嵌套，行列嵌套)
    int round = 0;
    int *rnd = (int *)malloc(sizeof(int) * NSIZE);
    memset(rnd, -1, sizeof(int) * NSIZE);

    // 两种处理方式：全部顶点：轮次。行顶点：轮次，列顶点：轮次。是否等价。

    auto start_kasi = std::chrono::system_clock::now();

    int x = 0;

    int *stack = (int *)malloc(sizeof(int) * NSIZE);
    int stack_last;

    // 主循环，直到处理完所有的边? 并且度数一和度数二的桶都为空
    while (!bucket1.empty() || !bucket2.empty() || !bucket3.empty())
    {
        // 度数一的顶点
        while (!bucket1.empty())
        {
            int v = bucket1.top();
            bucket1.pop();
            if (rg->remdegree[v] == 1)
            {
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
                basic_remove_neighborhood(u, rg->edge[u], rg->remdegree, bucket1, bucket2); // remdegree和实际的size对不上？？？   后添加的列没有到指定的行
                rg->edge[v].clear();
                deg1count++;
            }
        }

        //==================================================================================================================================================//
        //=======================================按照轮次，平衡处理====================================================================================//
        //==================================================================================================================================================//

        // 规则应用的有问题，存在可以应用规则的顶点没有应用
        if (!bucket2.empty())
        {
            int v = bucket2.top();
            bucket2.pop();

            x++;
            // if (x == 1639582)
            // {
            //     cout << endl;
            // }

            while (rnd[v] == round && !bucket2.empty()) // 直到找到一个本轮可以处理的顶点
            {
                bucket3.push(v);
                v = bucket2.top();
                bucket2.pop();
            }
            if (rg->remdegree[v] == 2)
            {
                //==================================================================================================================================================//
                //==============================================================搜索==================================================================================//
                //==================================================================================================================================================//
                row_que_start = 0, row_que_end = 0, col_que_start = 0, col_que_end = 0;

                row_que[row_que_end++] = v;
                visited[2][v] = v;
                // 向量中可能存储了度数为0的顶点，使用size()可能错误并且产生消耗大
                int point = 0;
                while (!rg->remdegree[rg->edge[v][point].endpoint])
                    point++;
                col_que[col_que_end++] = rg->edge[v][point].endpoint;
                visited[2][rg->edge[v][point].endpoint] = v;

                Edge e1 = rg->edge[v][point];
                point++;
                while (!rg->remdegree[rg->edge[v][point].endpoint]) // bucket2中的顶点的度数没有正常更新？？？？？？？？？？？
                    point++;
                col_que[col_que_end++] = rg->edge[v][point].endpoint;
                visited[2][rg->edge[v][point].endpoint] = v;
                Edge e2 = rg->edge[v][point];

                augment_recover_graph(e1, rg);
                augment_recover_graph(e2, rg);
                add_twins(e1, e2, rg); // 每次处理一个pair直接构建图

                while (col_que_start != col_que_end)
                {
                    int quecol = col_que[col_que_start++];

                    for (int ptr = 0; ptr < rg->edge[quecol].size(); ptr++)
                    {
                        int temp_row = rg->edge[quecol][ptr].endpoint;

                        if (rg->remdegree[temp_row] == 0 || visited[2][temp_row] == v || rnd[temp_row] == round) // 是否要判断度数                             访问的每一个行顶点要判断轮次
                        {
                            continue;
                        }
                        e1 = rg->edge[quecol][ptr]; // pair的第一条边
                        // 对于可以处理的顶点，更新其计数
                        if (visited[0][temp_row] != v)
                        {
                            visited[0][temp_row] = v;
                            visited[1][temp_row] = 1;
                        }
                        else
                        {
                            visited[1][temp_row]++;
                        }

                        // 判断计数是否满足条件加入队列
                        if (rg->remdegree[temp_row] - visited[1][temp_row] <= 1)
                        {
                            row_que[row_que_end++] = temp_row;
                            visited[2][temp_row] = v;

                            // 如果恰好满足要求，查找对应列
                            if (rg->remdegree[temp_row] - visited[1][temp_row] == 1)
                            {
                                for (int rowptr = 0; rowptr < rg->edge[temp_row].size(); rowptr++)
                                {
                                    int temp_col = rg->edge[temp_row][rowptr].endpoint;

                                    // 没有被删除的，没有加入过的
                                    if (rg->remdegree[temp_col] != 0)
                                    {
                                        if (visited[2][temp_col] != v) // 最后被加入的一条边
                                        {
                                            col_que[col_que_end++] = temp_col;
                                            visited[2][temp_col] = v;

                                            // 往前找一条边
                                            augment_recover_graph(e1, rg);
                                            e2 = rg->edge[temp_row][rowptr];
                                            augment_recover_graph(e2, rg);
                                            add_twins(e1, e2, rg); // 每次处理一个pair直接构建图
                                        }
                                    }
                                }
                            }
                            if (row_que_end == col_que_end) // 行列相等怎么处理
                            {
                                break;
                            }
                        }
                    }
                    if (row_que_end == col_que_end)
                    {
                        break;
                    }
                }

                //==================================================================================================================================================//
                //===================================================更新缩减图======================================================================================//
                //==================================================================================================================================================//

                // 删除行顶点
                if (row_que_end == col_que_end)
                {
                    row_que_end--;
                }

                deg2count += row_que_end;

                for (int i = 0; i < row_que_end; i++)
                {
                    rg->remdegree[row_que[i]] = 0;
                    rg->edge[row_que[i]].clear();
                }

                // 处理列顶点
                // 记录第一个列顶点

                // 寻找具有最大度数的列顶点
                int max_vtx_idx = 0;
                for (int i = 1; i < col_que_end; i++)
                {
                    if (rg->remdegree[col_que[i]] > rg->remdegree[col_que[max_vtx_idx]])
                    {
                        max_vtx_idx = i;
                    }
                }
                swap(col_que[0], col_que[max_vtx_idx]);

                // // 更新该列顶点的邻居
                // int upper_limit = rg->edge[col_que[0]].size() - 1; // 这种使用最后元素代替当前元素的方法只能够从后往前处理。    //从前往后搜素处理一次
                // for (int i = upper_limit; i >= 0; --i)
                // {
                //     if (rg->remdegree[rg->edge[col_que[0]][i].endpoint]) // 更新行顶点
                //     {
                //         mergearray[rg->edge[col_que[0]][i].endpoint] = col_que[0];
                //         rnd[rg->edge[col_que[0]][i].endpoint] = round;
                //     }
                //     else
                // remove_element(rg->edge[col_que[0]], i);
                // }

                int upper_limit = rg->edge[col_que[0]].size() - 1; // 这种使用最后元素代替当前元素的方法只能够从后往前处理。    //从前往后搜素处理一次
                stack_last = 0;
                for (int i = 0; i <= upper_limit; ++i)
                {
                    if (rg->remdegree[rg->edge[col_que[0]][i].endpoint]) // 更新行顶点
                    {
                        mergearray[rg->edge[col_que[0]][i].endpoint] = col_que[0];
                        rnd[rg->edge[col_que[0]][i].endpoint] = round;
                        stack[stack_last++] = i;
                    }
                }
                int temp_record = 0;
                for (int i = 0; i < stack_last; i++)
                {
                    rg->edge[col_que[0]][temp_record] = rg->edge[col_que[0]][stack[i]];
                    temp_record++;
                }
                for (int i = 0; i < upper_limit + 1 - stack_last; i++)
                {
                    rg->edge[col_que[0]].pop_back();
                }

                // 有的被移除的元素是此次被合并后才被移除的行顶点，使用--无法确定数量

                // 处理之后的每一个列顶点
                for (int i = 1; i < col_que_end; i++)
                {
                    // 两种处理方式，第一种处理完所有存在的边（通过度数记录）第二种处理完所有记录的边（直接访问到size）。
                    upper_limit = rg->edge[col_que[i]].size() - 1;
                    for (int j = 0; j <= upper_limit; ++j)
                    {
                        if (rg->remdegree[rg->edge[col_que[i]][j].endpoint])
                        {
                            rnd[rg->edge[col_que[i]][j].endpoint] = round;
                            if (mergearray[rg->edge[col_que[i]][j].endpoint] != col_que[0])
                            {
                                mergearray[rg->edge[col_que[i]][j].endpoint] = col_que[0];
                                // 从被合并顶点指向邻居顶点的端点就是邻居顶点，不需要修改，直接插入被合并顶点即可
                                rg->edge[col_que[0]].push_back(rg->edge[col_que[i]][j]);
                                // 连接cb边

                                // 连接bc边有点问题
                                int vtx_b = rg->edge[col_que[i]][j].endpoint;
                                for (int k = 0; k < rg->edge[vtx_b].size(); ++k)
                                {
                                    if (rg->edge[vtx_b][k].endpoint == col_que[i])
                                    {
                                        rg->edge[vtx_b][k].endpoint = col_que[0];
                                        rg->edge[vtx_b][k].original_start = rg->edge[col_que[i]][j].original_end;
                                        rg->edge[vtx_b][k].original_end = rg->edge[col_que[i]][j].original_start;
                                        break;
                                    }
                                }
                                // tmpedge.endpoint = col_que[0];
                                // tmpedge.original_start = rg->edge[col_que[i]][j].original_end;
                                // tmpedge.original_end = rg->edge[col_que[i]][j].original_start;
                                // rg->edge[rg->edge[col_que[i]][j].endpoint].push_back(tmpedge);
                            }
                            else
                            {
                                reduce_degree(rg->edge[col_que[i]][j].endpoint, rg->remdegree, bucket1, bucket3); // 子图中顶点可能被修改多次度数，直接加入bucket可能出错？直接加入，取出的时候进行一次判断
                            }
                        }
                    }

                    // 删除被合并列顶点
                    rg->remdegree[col_que[i]] = 0;
                    rg->edge[col_que[i]].clear();
                }

                // 按照轮次处理不处理这种顶点！！       如果一开始度数信息是正确的话，直接使用--是对的，这里使用了size，因此前面不需要更新度数信息
                rg->remdegree[col_que[0]] = rg->edge[col_que[0]].size();
                if (rg->remdegree[col_que[0]] == 1)
                {
                    bucket1.push(col_que[0]);
                }
                else if (rg->remdegree[col_que[0]] == 2)
                {
                    bucket3.push(col_que[0]);
                }
            }
        }
        else // 如果桶二为空，将桶三中元素全部置入bucket2
        {
            round++;
            swap(bucket2, bucket3);
        }
    }

    auto end_kasi = std::chrono::system_clock::now();
    cout << "存储的匹配基数为   " << deg1count + deg2count << " 其中deg1数量为 " << deg1count << " deg2数量为 " << deg2count << endl;

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // of << deg1count + deg2count << "    " << endl;
    // // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // // of << elapsed_seconds_kasi.count() << endl;
    // of.close();
    // 精确算法
    // 对于度数为1的顶点，直接设置相关的边的顶点的度数为0表示删除，并且更新相关顶点的度数没有实际更新边表。
    // 对于度数为2的顶点，修改其中一个顶点的连接状态，删除另外两个顶点。
    // 度数为0的顶点表示被删除

    free(mergearray);
    free(row_que);
    free(col_que);
    free(visited[0]);
    free(visited[1]);
    free(visited[2]);
    free(visited);
    free(stack);
}
