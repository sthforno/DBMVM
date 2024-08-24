#ifndef _DGMVM
#define _DGMVM
#include <cstring>
#include <fstream>

#include "kasi.h"
// 创新点一：结合间隙与连接的方法，保证插入效率的同时保证时间复杂度，是否可以单纯使用连接的方法？？？？？
// 创新点二：根据算法特点，限制间隙的大小从而保证空间复杂度
// 创新点三：单次插入和批量插入优化
// 潜在创新点：每次插边应该插入两个边，现有基于连接的方法没有考察插入和遍历频繁交替的情况。

// 动态图存储结构    间隙加上链接的方法

// 问题：一个度数为一的顶点v与一个度数很高的顶点u相连，u一开始的度数不高，后续不断增加。

struct dynamic_graph
{
    int n, m, edge_size, times;
    int *vtx_ptr_start;
    int *vtx_ptr_end;
    int *cur_edge;          // 当前指向
    int *origin_edge_start; // 原始的指向
    int *origin_edge_end;   // 原始的终点
    int *matching;
    int *degree;
    int *vtx_link_next; // 存储下一个边表
    int *vtx_link_last; // 存储最后一个边表
    int *vtx_link_temp; // 存储第一个空表
};

// 匹配树存储结构
struct match_tree
{
    int *project;
    int *re_project;

    // 记录所有在树上的顶点
    int *tree_vtx;
    int tree_vtx_size;

    // 记录树上顶点的度数
    int *tree_vtx_degree;

    int *tree_vtx_ptr_start;
    int *tree_vtx_ptr_end;
    int *tree_edge;
};
void DGMVM(graph *og, dynamic_graph *dg, match_tree *mt, int times);
void DBMVM(graph *og, dynamic_graph *dg, match_tree *mt, int times);
inline void init_dynamic_graph(graph *og, dynamic_graph *dg, stack<int> &bucket1, stack<int> &bucket2, int times) // 如果在h文件中实现，一定要加上inline
{
    dg->n = og->n;
    dg->m = og->m;
    dg->times = times;
    dg->edge_size = dg->times * dg->m;
    dg->vtx_ptr_start = (int *)malloc(sizeof(int) * (dg->n + 1)); // 使用两个数组存储顶点的边索引位置比较好还是一个比较好？？？？？？
    dg->vtx_ptr_end = (int *)malloc(sizeof(int) * dg->n);
    dg->cur_edge = (int *)malloc(sizeof(int) * dg->edge_size);
    dg->origin_edge_start = (int *)malloc(sizeof(int) * dg->edge_size);
    dg->origin_edge_end = (int *)malloc(sizeof(int) * dg->edge_size);
    memset(dg->cur_edge, -1, sizeof(int) * dg->edge_size);
    memset(dg->origin_edge_start, -1, sizeof(int) * dg->edge_size);
    memset(dg->origin_edge_end, -1, sizeof(int) * dg->edge_size);

    dg->degree = (int *)malloc(sizeof(int) * dg->n);
    dg->vtx_link_next = (int *)malloc(sizeof(int) * dg->n);
    memset(dg->vtx_link_next, -1, sizeof(int) * dg->n);
    dg->vtx_link_last = (int *)malloc(sizeof(int) * dg->n);
    dg->vtx_link_temp = (int *)malloc(sizeof(int) * dg->n);

    for (int i = 0; i < dg->n; i++) // 复制图   优化点一：给顶点分配边数大小具有上界n
    {
        dg->vtx_ptr_start[i] = dg->times * og->vtx_pointer[i];
        dg->vtx_ptr_end[i] = dg->vtx_ptr_start[i];
        dg->vtx_link_temp[i] = i;
        dg->vtx_link_last[i] = i;
        for (int ptr = og->vtx_pointer[i]; ptr < og->vtx_pointer[i + 1]; ptr++)
        {
            dg->cur_edge[dg->vtx_ptr_end[i]] = og->endV[ptr];
            dg->origin_edge_end[dg->vtx_ptr_end[i]] = og->endV[ptr];
            dg->origin_edge_start[dg->vtx_ptr_end[i]] = i;
            dg->vtx_ptr_end[i]++;
        }
        dg->degree[i] = og->vtx_pointer[i + 1] - og->vtx_pointer[i];
        if (dg->degree[i] == 1)
        {
            bucket1.push(i);
        }
        else if (dg->degree[i] == 2)
        {
            bucket2.push(i);
        }
    }

    dg->vtx_ptr_start[dg->n] = dg->edge_size;

    dg->matching = (int *)malloc(sizeof(int) * og->n);
    memset(dg->matching, -1, sizeof(int) * dg->n);
}

inline void free_dynamic_graph(dynamic_graph *dg)
{
    free(dg->vtx_ptr_start);
    free(dg->vtx_ptr_end);
    free(dg->cur_edge);
    free(dg->origin_edge_end);
    free(dg->origin_edge_start);
    free(dg->degree);
    free(dg->matching);
    free(dg->vtx_link_next);
    free(dg->vtx_link_last);
    free(dg->vtx_link_temp);
    // free(dg);
}

// 删除顶点是否需要更新顶点连接的所有的边表的状态？？？？？？？？？？？？？？？？？？？
inline void delet_vertex(int v, dynamic_graph *dg)
{
    dg->degree[v] = 0;
    dg->vtx_ptr_end[v] = dg->vtx_ptr_start[v];
    int iter_vtx = v;
    int temp_vtx;
    while (dg->vtx_link_next[iter_vtx] != -1)
    {
        temp_vtx = dg->vtx_link_next[iter_vtx];                  // 获取新的顶点
        dg->vtx_link_next[iter_vtx] = -1;                        // 将当前迭代顶点移除
        iter_vtx = temp_vtx;                                     // 更新迭代顶点
        dg->vtx_ptr_end[iter_vtx] = dg->vtx_ptr_start[iter_vtx]; // 重置新的迭代顶点的访问范围
    }
}
// 只有在插入边的时候可能需要平衡整个数据结构,
// 每次插入边要在两个位置插入边
// 四种情况 1）ab都存在空隙 2）a存在空隙 3)b存在空隙 4）ab都不存在空隙。
// 插边的优化策略是否需要同时考虑这两个位置？？？？？？？？？？？？？？？？？？？？？？？？

// 将ab边改为cb边
// 传入ab边的edge_id，传入a，b，c三个顶点
// bug 1 插边没有更新度数
// bug 2 第一次连接到新的表时，该表应该是空的，直接往里存东西的话，会记录重复的边
// bug 3 错误使用变量名，导致没有正确更新边表的末尾
// bug 4 优先处理度数为一的顶点，将while改为if
// bug 5 没有确定边界条件，边表最后一个位置是否以及连接了
// bug 6 错误使用变量，导致移除了特定的边
// bug 7 混淆了指针，边表最后一个元素存储的是顶点而不是顶点指向的边表
// bug 8 边界设置有问题，表示连接到边表，该值应该大于最大的边表元素，dg->m  ---> 2*dg->m
// bug 9 存在长边表被合并移除的情况，此时它自己无法存储他自己。
// bug 10 直接连接的话，边表中可能会记录重复的元素
// bug 11 将边表中当前元素与最后一个元素互换之后，没有处理最后一个元素
// bug 12 把当前顶点需要记录的信息先记录下来再移除
// bug 13 表现，将间隙由大于等于改为大于后，添加了连接的概率，增加匹配基数。 ?????????
// bug 14 边界条件问题，边表可能只存在一个元素，移除该元素之后，边表指针指到前面一个边表去了
// bug 15 改变遍历方法之后，之前的跳出循环的方法可能不能够再跳出循环
// bug 16 使用memcpy出错的原因，变量的类型不一样，比如int和long

// 将vtx_a与顶点栈中的连边全部转换为vtx_c与其相连
inline void insert_edgelist_to_dynamic_graph(int vtx_a, int *stack_vtx, int *stack_edge_id, int stack_end, int vtx_c, dynamic_graph *dg)
{
    // 移除顶点a
    delet_vertex(vtx_a, dg);
    // 连接bc边
    int iter_vtx;
    // 处理每一个行顶点
    for (int i = 0; i < stack_end; i++)
    {
        iter_vtx = stack_vtx[i]; // 寻找到一个有空隙的边表
        int curret_vtx = dg->vtx_link_temp[iter_vtx];
        while (dg->vtx_ptr_end[curret_vtx] == dg->vtx_ptr_start[curret_vtx + 1] && curret_vtx != -1)
        {
            curret_vtx = dg->vtx_link_next[curret_vtx];
        }
        dg->vtx_link_temp[iter_vtx] = curret_vtx;

        if (curret_vtx != -1) // 有可能找到的顶点是没有空隙的
        {
            dg->cur_edge[dg->vtx_ptr_end[curret_vtx]] = vtx_c;
            dg->origin_edge_start[dg->vtx_ptr_end[curret_vtx]] = dg->origin_edge_start[stack_edge_id[i]];
            dg->origin_edge_end[dg->vtx_ptr_end[curret_vtx]] = dg->origin_edge_end[stack_edge_id[i]];
            dg->vtx_ptr_end[curret_vtx]++;

            if (dg->vtx_ptr_end[curret_vtx] == dg->vtx_ptr_start[curret_vtx + 1]) // 插满了将间隙指针往后移
            {
                dg->vtx_link_temp[stack_vtx[i]] = dg->vtx_link_next[curret_vtx];
            }
        }
        else // 找不到直接更新所有边表
        {
            do
            {
                for (int ptr = dg->vtx_ptr_start[iter_vtx]; ptr < dg->vtx_ptr_end[iter_vtx]; ptr++)
                {
                    int temp_vtx = dg->cur_edge[ptr];
                    if (dg->degree[temp_vtx] == 0) // 如果顶点不存在
                    {
                        int last_ptr = dg->vtx_ptr_end[iter_vtx] - 1; // 还可以继续优化,直接找到存在的顶点！！！！！！！！！！！！！！！！

                        dg->cur_edge[ptr] = dg->cur_edge[last_ptr];
                        dg->origin_edge_start[ptr] = dg->origin_edge_start[last_ptr];
                        dg->origin_edge_end[ptr] = dg->origin_edge_end[last_ptr];

                        // 更新末尾并且重新处理当前顶点
                        dg->vtx_ptr_end[iter_vtx]--;
                        if (dg->vtx_ptr_end[iter_vtx] < dg->vtx_ptr_start[iter_vtx])
                        {
                            cout << "wrong" << endl;
                        }
                        ptr--;
                    }
                }
                iter_vtx = dg->vtx_link_next[iter_vtx];
            } while (iter_vtx != -1);

            iter_vtx = stack_vtx[i];
            while (dg->vtx_ptr_end[iter_vtx] == dg->vtx_ptr_start[iter_vtx + 1] && iter_vtx != -1)
            {
                iter_vtx = dg->vtx_link_next[iter_vtx];
            }

            if (iter_vtx != -1)
            {
                dg->cur_edge[dg->vtx_ptr_end[iter_vtx]] = vtx_c;
                dg->origin_edge_start[dg->vtx_ptr_end[iter_vtx]] = dg->origin_edge_start[stack_edge_id[i]];
                dg->origin_edge_end[dg->vtx_ptr_end[iter_vtx]] = dg->origin_edge_end[stack_edge_id[i]];
                dg->vtx_ptr_end[iter_vtx]++;

                if (dg->vtx_ptr_end[iter_vtx] == dg->vtx_ptr_start[iter_vtx + 1]) // 插满了将间隙指针往后移
                {
                    dg->vtx_link_temp[stack_vtx[i]] = dg->vtx_link_next[iter_vtx];
                }
            }
            else
            {
                cout << "wrong" << endl;
            }

            // 更新边表之后再进行插入
        }
    }
    // 连接cb边

    iter_vtx = dg->vtx_link_last[vtx_c];

    if (iter_vtx == -1)
    {
        iter_vtx = vtx_c;
    }

    for (int i = 0; i < stack_end; i++) // 将栈中的顶点全部连接到b上，并且更新顶点b的度数
    {
        dg->cur_edge[dg->vtx_ptr_end[iter_vtx]] = stack_vtx[i];
        dg->origin_edge_start[dg->vtx_ptr_end[iter_vtx]] = dg->origin_edge_start[stack_edge_id[i]];
        dg->origin_edge_end[dg->vtx_ptr_end[iter_vtx]] = dg->origin_edge_end[stack_edge_id[i]];
        dg->vtx_ptr_end[iter_vtx]++;
    }
    dg->degree[vtx_c] += stack_end;
}

// 间隙不够的情况，直接连接，使用连接的方法，移除a顶点只重置度数
inline void connect_edgelist_to_dynamic_graph(int vtx_a, int *stack_vtx, int *stack_edge_id, int stack_end, int vtx_c, dynamic_graph *dg)
{
    // 删除顶点a
    dg->degree[vtx_a] = 0;
    // 连接bc边
    int iter_vtx;
    // 处理每一个行顶点
    for (int i = 0; i < stack_end; i++)
    {
        iter_vtx = stack_vtx[i]; // 寻找到一个有空隙的顶点
        int curret_vtx = dg->vtx_link_temp[iter_vtx];
        while (dg->vtx_ptr_end[curret_vtx] == dg->vtx_ptr_start[curret_vtx + 1] && curret_vtx != -1)
        {
            curret_vtx = dg->vtx_link_next[curret_vtx];
        }
        dg->vtx_link_temp[iter_vtx] = curret_vtx;

        if (curret_vtx != -1) // 找到空隙直接插
        {
            dg->cur_edge[dg->vtx_ptr_end[curret_vtx]] = vtx_c;
            dg->origin_edge_start[dg->vtx_ptr_end[curret_vtx]] = dg->origin_edge_start[stack_edge_id[i]];
            dg->origin_edge_end[dg->vtx_ptr_end[curret_vtx]] = dg->origin_edge_end[stack_edge_id[i]];
            dg->vtx_ptr_end[curret_vtx]++;

            if (dg->vtx_ptr_end[curret_vtx] == dg->vtx_ptr_start[curret_vtx + 1]) // 插满了将间隙指针往后移
            {
                dg->vtx_link_temp[stack_vtx[i]] = dg->vtx_link_next[curret_vtx];
            }
        }
        else // 找不到直接更新所有边表
        {
            do
            {
                for (int ptr = dg->vtx_ptr_start[iter_vtx]; ptr < dg->vtx_ptr_end[iter_vtx]; ptr++)
                {
                    int temp_vtx = dg->cur_edge[ptr];
                    if (dg->degree[temp_vtx] == 0 || temp_vtx == -1) // 如果顶点不存在
                    {
                        int last_ptr = dg->vtx_ptr_end[iter_vtx] - 1; // 还可以继续优化,直接找到存在的顶点！！！！！！！！！！！！！！！！

                        dg->cur_edge[ptr] = dg->cur_edge[last_ptr];
                        dg->origin_edge_start[ptr] = dg->origin_edge_start[last_ptr];
                        dg->origin_edge_end[ptr] = dg->origin_edge_end[last_ptr];

                        // 更新末尾并且重新处理当前顶点
                        dg->vtx_ptr_end[iter_vtx]--;
                        ptr--;
                    }
                }
                iter_vtx = dg->vtx_link_next[iter_vtx];
            } while (iter_vtx != -1);

            // 在更新后的边表中找到一个间隙
            iter_vtx = stack_vtx[i];
            while (dg->vtx_ptr_end[iter_vtx] == dg->vtx_ptr_start[iter_vtx + 1] && iter_vtx != -1)
            {
                iter_vtx = dg->vtx_link_next[iter_vtx];
            }

            if (iter_vtx != -1)
            {
                dg->cur_edge[dg->vtx_ptr_end[iter_vtx]] = vtx_c;
                dg->origin_edge_start[dg->vtx_ptr_end[iter_vtx]] = dg->origin_edge_start[stack_edge_id[i]];
                dg->origin_edge_end[dg->vtx_ptr_end[iter_vtx]] = dg->origin_edge_end[stack_edge_id[i]];
                dg->vtx_ptr_end[iter_vtx]++;

                if (dg->vtx_ptr_end[iter_vtx] == dg->vtx_ptr_start[iter_vtx + 1]) // 插满了将间隙指针往后移
                {
                    dg->vtx_link_temp[stack_vtx[i]] = dg->vtx_link_next[iter_vtx];
                }
            }
            else
            {
                cout << "wrong" << endl;
            }

            // 更新边表之后再进行插入
        }
    }
    // 连接cb边
    iter_vtx = dg->vtx_link_last[vtx_c];
    dg->vtx_link_next[iter_vtx] = vtx_a;                 // 连接最后一个边表
    dg->vtx_link_last[vtx_c] = dg->vtx_link_last[vtx_a]; // 更新指向的最后一个边表

    dg->degree[vtx_c] += stack_end;
}

inline void init_match_tree(match_tree *mt, graph *dg)
{
    mt->project = (int *)malloc(sizeof(int) * dg->n);
    memset(mt->project, -1, sizeof(int) * dg->n);
    mt->re_project = (int *)malloc(sizeof(int) * dg->n);
    memset(mt->re_project, -1, sizeof(int) * dg->n);

    mt->tree_vtx = (int *)malloc(sizeof(int) * dg->n);
    mt->tree_vtx_size = 0;

    // 存储树上的边
    mt->tree_vtx_ptr_start = (int *)malloc(sizeof(int) * (dg->n + 1));
    mt->tree_vtx_ptr_end = (int *)malloc(sizeof(int) * (dg->n + 1));
    for (int i = 0; i <= dg->n; i++)
    {
        mt->tree_vtx_ptr_start[i] = dg->vtx_pointer[i];
        mt->tree_vtx_ptr_end[i] = dg->vtx_pointer[i];
    }
    mt->tree_edge = (int *)malloc(sizeof(int) * dg->m);

    mt->tree_vtx_degree = (int *)malloc(sizeof(int) * dg->n);
    memset(mt->tree_vtx_degree, 0, sizeof(int) * dg->n);
}

inline void free_match_tree(match_tree *mt)
{
    free(mt->project);
    free(mt->re_project);
    free(mt->tree_edge);
    free(mt->tree_vtx);
    free(mt->tree_vtx_degree);
    free(mt->tree_vtx_ptr_start);
    free(mt->tree_vtx_ptr_end);
}

// 将被合并的边记录到树上，包含顶点，顶点度数以及顶点之间的边连接
inline void augment_recover_graph(int insert_edge_id, match_tree *mt, dynamic_graph *dg)
{
    int original_a = dg->origin_edge_start[insert_edge_id];
    int original_b = dg->origin_edge_end[insert_edge_id];

    // 记录顶点       记录边的端点并且更新其在树上的度数
    if (mt->tree_vtx_degree[original_a] == 0) // this is done to save some space
        mt->tree_vtx[mt->tree_vtx_size++] = original_a;
    mt->tree_vtx_degree[original_a]++;
    if (mt->tree_vtx_degree[original_b] == 0)
        mt->tree_vtx[mt->tree_vtx_size++] = original_b;
    mt->tree_vtx_degree[original_b]++;

    // 记录边        恢复图的存储仍然使用CSR类似的格式
    mt->tree_edge[mt->tree_vtx_ptr_end[original_b]++] = original_a;
    mt->tree_edge[mt->tree_vtx_ptr_end[original_a]++] = original_b;
}

// 恢复匹配
inline void recover_matching_from_match_tree(graph *og, dynamic_graph *dg, match_tree *mt, long *mate)
{
    // 将缩减图的匹配转变为全局匹配
    // int temp_num1 = 0, temp_num2 = 0, temp_num3 = 0, temp_edge_num = 0;
    int iter_vtx;
    bool find;
    for (int i = 0; i < og->n; i++)
    {
        if (mate[i] != -1)
        {
            iter_vtx = mt->re_project[i];
            find = false;
            do
            {
                for (int ptr = dg->vtx_ptr_start[iter_vtx]; ptr < dg->vtx_ptr_end[iter_vtx]; ptr++)
                {
                    if (dg->cur_edge[ptr] == mt->re_project[mate[i]])
                    {
                        find = true;
                        dg->matching[dg->origin_edge_start[ptr]] = dg->origin_edge_end[ptr];
                        dg->matching[dg->origin_edge_end[ptr]] = dg->origin_edge_start[ptr];
                        break;
                    }
                }
                iter_vtx = dg->vtx_link_next[iter_vtx];
            } while (iter_vtx != -1 && find == false);

            // dg->matching[mt->re_project[i]] = mt->re_project[mate[i]];
            // dg->matching[mt->re_project[mate[i]]] = mt->re_project[i];
            // temp_num1++;
        }
    }

    // int incr = 0;
    // 将存储在匹配树上的顶点分为两类
    stack<int> bucket1_recover;
    stack<int> bucket_matched_recover;
    for (int k = 0; k < mt->tree_vtx_size; ++k)
    {
        // temp_edge_num += mt->tree_vtx_degree[mt->tree_vtx[k]];
        if (mt->tree_vtx_degree[mt->tree_vtx[k]] == 1)
        {
            bucket1_recover.push(mt->tree_vtx[k]);
            // temp_num2++;
        }
        if (dg->matching[mt->tree_vtx[k]] != -1)
        {
            bucket_matched_recover.push(mt->tree_vtx[k]);
            // temp_num3++;
        }
    }

    // cout << temp_num1 << " " << temp_num2 << " " << temp_num3 << " " << temp_edge_num << endl;

    // int card_of_recover_matching = 0;

    // int card_of_matching1 = 0;
    // for (int i = 0; i < dg->n; i++)
    // {
    //     if (dg->matching[i] != -1)
    //     {
    //         card_of_matching1++;
    //     }
    // }
    // cout << "恢复之前的匹配为 " << card_of_matching1 << endl;

    // 将这两类顶点处理，直到为空，先处理与外界相连的顶点
    while (!bucket1_recover.empty() || !bucket_matched_recover.empty())
    {
        while (!bucket_matched_recover.empty())
        {
            // 取出被匹配的树顶点以及其邻居
            int v = bucket_matched_recover.top();
            int mvv = dg->matching[v];
            bucket_matched_recover.pop();

            // 处理被匹配顶点的每一个邻居
            for (int ptr = mt->tree_vtx_ptr_start[v]; ptr < mt->tree_vtx_ptr_end[v]; ++ptr)
            {
                int u = mt->tree_edge[ptr];
                if (u != mvv)
                { // 如果邻居不是当前顶点的匹配顶点，并且存在，更新其度数
                    if (mt->tree_vtx_degree[u] > 0)
                    {
                        mt->tree_vtx_degree[u]--;
                        if (mt->tree_vtx_degree[u] == 1)
                            bucket1_recover.push(u);
                    }
                }
            }
            mt->tree_vtx_degree[v] = 0;
        }
        // 处理匹配树中的叶子顶点
        if (!bucket1_recover.empty())
        {
            int v = bucket1_recover.top();
            bucket1_recover.pop();

            // 如果叶子顶点还存在并且没有被处理过
            if (mt->tree_vtx_degree[v] > 0 && dg->matching[v] == -1)
            {

                // 寻找满足要求的邻居顶点
                int temp_neigh;
                for (int ptr = mt->tree_vtx_ptr_start[v]; ptr < mt->tree_vtx_ptr_end[v]; ptr++)
                {
                    temp_neigh = mt->tree_edge[ptr];
                    // 没有匹配，并且存在，并且满足当前匹配边的对应边没有被匹配（？？？？？？？？？？？？？？？？？？？？？？，是否需要这个条件）
                    if (dg->matching[temp_neigh] == -1 && mt->tree_vtx_degree[temp_neigh] != 0) //&& dg->matching[mt->treedge_twin_1[v][neigh]] != mt->treedge_twin_2[v][neigh]
                    {
                        break;
                    }
                }

                // 将匹配的两个顶点加入被匹配的叶子顶点集合
                dg->matching[v] = temp_neigh;
                dg->matching[temp_neigh] = v;
                // incr++;
                bucket_matched_recover.push(temp_neigh);
                bucket_matched_recover.push(v);

                mt->tree_vtx_degree[v] = 0;

                // card_of_recover_matching++;
            }
        }
    }

    int card_of_matching = 0;
    for (int i = 0; i < dg->n; i++)
    {
        if (dg->matching[i] != -1)
        {
            card_of_matching++;
        }
    }
    cout << "The cardinality of the reconstructed matching: " << card_of_matching << endl;

    // cout << "从存储匹配中恢复的匹配数量为" << card_of_recover_matching << endl;

    // std::ofstream of;
    // of.open("TIME.txt", ios ::app);
    // of << card_of_matching << "  ";
    // of.close();

    // free_match_tree(mt);

    free(mate);
}

// 根据度数信息缩减og
inline void reduce_the_dynamic_graph(graph *og, dynamic_graph *dg, match_tree *mt)
{
    // 修改og,首先记录存在的顶点的数量
    int num_of_vtx = 0;
    for (int i = 0; i < dg->n; ++i)
    {
        if (dg->degree[i] != 0)
        {
            mt->project[i] = num_of_vtx;
            mt->re_project[num_of_vtx] = i;
            num_of_vtx++;
        }
    }

    // 获取行顶点的数量
    int nrows = 0;
    for (int i = og->nrows; i < dg->n; ++i)
    {
        if (mt->project[i] != -1)
        {
            nrows = mt->project[i];
            break;
        }
    }
    og->nrows = nrows;

    // 重新构建og
    // 处理每一个非零元素
    int iter_vtx;
    int neighbor_vtx;
    int num_of_edge;
    int examined;
    int j;

    // 将每一个顶点的边记录到新的og中去
    for (int i = 0; i < num_of_vtx; i++)
    {
        iter_vtx = mt->re_project[i];     // 原始的顶点id
        num_of_edge = og->vtx_pointer[i]; // 边数
        do
        {
            for (int ptr = dg->vtx_ptr_start[iter_vtx]; ptr < dg->vtx_ptr_end[iter_vtx]; ptr++)
            {
                if (dg->degree[dg->cur_edge[ptr]])
                {
                    og->endV[num_of_edge++] = mt->project[dg->cur_edge[ptr]];
                }
            }
            iter_vtx = dg->vtx_link_next[iter_vtx];
        } while (iter_vtx != -1);
        og->vtx_pointer[i + 1] = num_of_edge;
    }

    // cout << "总边数   " << num_of_edge << "顶点数   " << num_of_vtx << endl;

    og->n = num_of_vtx;
    og->m = num_of_edge;
}

// 移除顶点后要重置end指针以及该边表中的所有元素!!!!!!!!!!!!!!!!

// 合并边操作分为以下几个步骤
// 1）插入边 将被移除邻居顶点与其他顶点相关的边进行插入
// 1.2  将ab边    变为cb边       在c顶点的边表中插入ab边     b顶点的边表中将原始的ab变为cb边。
// 2）移除边  相关顶点与被移除的邻居顶点之间的边要移除
// 3）移除顶点 移除可合并顶点以及其中一个邻居顶点
// 4)

// 数据结构组织区别
// 使用两个数组记录顶点索引和使用一个数组记录顶点索引 ？？ 边的起始位置和终止位置
// 使用结构体插入和使用三个数组插入的效率？ 边的终止顶点和原始起点和原始终点

// 为什么要记录边的原始起点和终点？？？？？ 构建匹配树的时候需要，匹配的时候要恢复原始的匹配
#endif