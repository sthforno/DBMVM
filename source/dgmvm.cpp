#include <queue>
#include <stack>
#include <ctype.h>
#include <chrono>
#include <ctime>
#include <unordered_set>
#include <set>
#include <fstream>
#include <cstring>
#include "../include/kasi.h" //cpp文件里面可以随意包含
#include "../include/dgmvm.h"

// 遍历方法，while(iter<dgree)+for循环，for循环后统一更新iter。

// 测验是否有相同的边
void verify_vtx(int v, dynamic_graph *dg, int *stack, int *verify_array)
{

    int stack_last = 0;
    int iter_vtx = v;

    do
    {
        for (int ptr = dg->vtx_ptr_start[iter_vtx]; ptr < dg->vtx_ptr_end[iter_vtx]; ptr++)
        {
            if (dg->degree[dg->cur_edge[ptr]] && dg->cur_edge[ptr] != -1)
            {
                stack[stack_last++] = dg->cur_edge[ptr];
            }
        }
        iter_vtx = dg->vtx_link_next[iter_vtx];
    } while (iter_vtx != -1);

    //

    for (int i = 0; i < stack_last; i++)
    {
        if (verify_array[stack[i]] != v)
        {
            verify_array[stack[i]] = v;
        }
        else
        {
            for (int j = 0; j < i; j++)
            {
                if (stack[j] == stack[i])
                {
                    cout << "wrong" << endl;
                }
            }
            // 验证该边是否存在
        }
    }
}

// 测试两个顶点是否存在连接

void verify_edge_connection(int vtx_a, int vtx_b, dynamic_graph *dg)
{
    bool find1 = false;
    bool find2 = false;
    int iter_vtx = vtx_a;
    do
    {
        for (int ptr = dg->vtx_ptr_start[iter_vtx]; ptr < dg->vtx_ptr_end[iter_vtx]; ptr++)
        {
            if (dg->cur_edge[ptr] == vtx_b)
            {
                find1 = true;
                break;
            }
        }
        iter_vtx = dg->vtx_link_next[iter_vtx];
    } while (iter_vtx != -1);

    iter_vtx = vtx_b;
    do
    {
        for (int ptr = dg->vtx_ptr_start[iter_vtx]; ptr < dg->vtx_ptr_end[iter_vtx]; ptr++)
        {
            if (dg->cur_edge[ptr] == vtx_a)
            {
                find2 = true;
                break;
            }
        }
        iter_vtx = dg->vtx_link_next[iter_vtx];
    } while (iter_vtx != -1);

    if ((find1 == true && find2 == false) || (find2 == true && find1 == false))
    {
        cout << "出错了" << endl;
    }
}

// 显式的移除ab边消耗是大的，因为要从a表和b表中找到对方。
void DGMVM(graph *og, dynamic_graph *dg, match_tree *mt, int times)
{
    int u, w;
    int start_ptr, end_ptr;
    int start_ptr1, end_ptr1;
    int temp_vtx;
    int temp_vtx_neigh;

    int original_vtx_a, original_vtx_b;

    int iter_vtx;
    stack<int> bucket1;
    stack<int> bucket2;

    // rg->treeedge.resize(og->n);
    // rg->treedge_twin_1.resize(og->n);
    // rg->treedge_twin_2.resize(og->n);
    // rg->treevertexdegree = (int *)calloc(og->n, sizeof(int));

    // 初始化匹配

    int *mergearray = (int *)malloc(sizeof(int) * og->n);
    memset(mergearray, -1, sizeof(int) * og->n);

    // int *verify_array = (int *)malloc(sizeof(int) * og->n);
    // memset(verify_array, -1, sizeof(int) * og->n);

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // auto start_kasi = std::chrono::system_clock::now();
    init_dynamic_graph(og, dg, bucket1, bucket2, times);
    init_match_tree(mt, og);
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
    visited[0] = (int *)malloc(sizeof(int) * og->n);
    visited[1] = (int *)malloc(sizeof(int) * og->n);
    visited[2] = (int *)malloc(sizeof(int) * og->n);
    memset(visited[0], -1, sizeof(int) * og->n);
    memset(visited[1], -1, sizeof(int) * og->n);
    memset(visited[2], -1, sizeof(int) * og->n);

    // 存储需要移除的列顶点的邻居
    int *stack_vtx_neigh = (int *)malloc(sizeof(int) * og->n);
    int *stack_edge_id = (int *)malloc(sizeof(int) * og->n);
    int stack_last;

    int temp_num = 0;
    int deg1count = 0;
    int deg2count = 0;

    // 两种处理方式：全部顶点：轮次。行顶点：轮次，列顶点：轮次。是否等价。

    auto start_kasi = std::chrono::system_clock::now();

    int merge_operations = 0;
    // 主循环，直到处理完所有的边? 并且度数一和度数二的桶都为空
    while (!bucket1.empty() || !bucket2.empty()) //
    {
        // 度数一的顶点
        while (!bucket1.empty())
        {
            int v = bucket1.top();
            bucket1.pop();
            if (dg->degree[v] == 1) // 隐式的删除方法导致查找邻居需要访问邻接表中所有的元素
            {
                iter_vtx = v;
                u = -1;
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                    {
                        temp_vtx = dg->cur_edge[start_ptr];
                        if (dg->degree[temp_vtx] && temp_vtx != -1) // 是否需要大于0!!!!!!!!!!!!!!!!!!!!
                        {
                            u = temp_vtx;
                            break;
                        }
                    }

                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (u == -1 && iter_vtx != -1); // 找到u直接跳出循环

                // 将该边的原始边加入匹配
                dg->matching[dg->origin_edge_start[start_ptr]] = dg->origin_edge_end[start_ptr];
                dg->matching[dg->origin_edge_end[start_ptr]] = dg->origin_edge_start[start_ptr];

                // 移除顶点v
                delet_vertex(v, dg);

                // 更新顶点u相关的顶点

                iter_vtx = u; // 处理u的所有邻居
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++) // 访问当前顶点的边表
                    {
                        temp_vtx = dg->cur_edge[start_ptr];
                        // 更新tempvtx的状态，如何将tempvtx的邻居边表中u移除？ 只能隐式的移除，否则会增加时间复杂度
                        if (dg->degree[temp_vtx] && temp_vtx != -1) // 度数是准的，顶点指针不一定准
                        {
                            dg->degree[temp_vtx]--;
                            if (dg->degree[temp_vtx] == 1)
                                bucket1.push(temp_vtx);
                            else if (dg->degree[temp_vtx] == 2)
                                bucket2.push(temp_vtx);
                        }
                    }

                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (iter_vtx != -1);

                // 移除顶点u
                delet_vertex(u, dg);

                temp_num++;
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

            if (dg->degree[v] == 2)
            {
                merge_operations++;
                //==================================================================================================================================================//
                //==============================================================搜索==================================================================================//
                //==================================================================================================================================================//
                row_que_start = 0, row_que_end = 0, col_que_start = 0, col_que_end = 0;

                row_que[row_que_end++] = v;
                visited[2][v] = v;

                // 查找顶点v的两个邻居顶点
                int edge_id[2];
                iter_vtx = v;
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                    {
                        temp_vtx = dg->cur_edge[start_ptr];
                        if (dg->degree[temp_vtx] && temp_vtx != -1)
                        {
                            edge_id[col_que_end] = start_ptr;
                            col_que[col_que_end++] = temp_vtx;
                            visited[2][temp_vtx] = v;
                        }
                    }
                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (col_que_end < 2 && iter_vtx != -1);

                augment_recover_graph(edge_id[0], mt, dg);
                augment_recover_graph(edge_id[1], mt, dg);
                // add_twins(e1, e2, rg); // 每次处理一个pair直接构建图

                while (col_que_start != col_que_end)
                {
                    int quecol = col_que[col_que_start++];

                    // 遍历quecol的每一个邻居
                    iter_vtx = quecol;
                    do
                    {
                        for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                        {
                            int temp_row = dg->cur_edge[start_ptr];

                            if (dg->degree[temp_row] == 0 || visited[2][temp_row] == v || temp_row == -1) // 是否要判断度数                             访问的每一个行顶点要判断轮次
                            {
                                continue;
                            }
                            edge_id[0] = start_ptr; // pair的第一条边
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
                            if (dg->degree[temp_row] - visited[1][temp_row] <= 1)
                            {
                                row_que[row_que_end++] = temp_row;
                                visited[2][temp_row] = v;

                                if (dg->degree[temp_row] - visited[1][temp_row] == 1)
                                {
                                    // 查找满足要求的行顶点的邻居顶点集合中没有被加入的列顶点
                                    int iter_vtx1 = temp_row;
                                    int find = false;
                                    do
                                    {
                                        for (start_ptr1 = dg->vtx_ptr_start[iter_vtx1]; start_ptr1 < dg->vtx_ptr_end[iter_vtx1]; start_ptr1++)
                                        {
                                            int temp_col = dg->cur_edge[start_ptr1];

                                            // 没有被删除的，没有加入过的
                                            if (dg->degree[temp_col] && temp_col != -1)
                                            {
                                                if (visited[2][temp_col] != v) // 最后被加入的一条边
                                                {
                                                    col_que[col_que_end++] = temp_col;
                                                    visited[2][temp_col] = v;
                                                    find = true;

                                                    // 插入匹配边需要信息
                                                    augment_recover_graph(start_ptr, mt, dg);
                                                    augment_recover_graph(start_ptr1, mt, dg);
                                                    // add_twins(e1, e2, rg);
                                                }
                                            }
                                        }

                                        iter_vtx1 = dg->vtx_link_next[iter_vtx1];
                                    } while (find == false && iter_vtx1 != -1);
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

                        iter_vtx = dg->vtx_link_next[iter_vtx];
                    } while (iter_vtx != -1);

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
                    delet_vertex(row_que[i], dg);
                }

                // 处理列顶点
                // 将度数最高的列顶点放在首位
                int max_vtx_idx = 0;
                for (int i = 1; i < col_que_end; i++)
                {
                    if (dg->degree[col_que[i]] > dg->degree[col_que[max_vtx_idx]])
                    {
                        max_vtx_idx = i;
                    }
                }
                swap(col_que[0], col_que[max_vtx_idx]);

                // 遍历具有最高度数的列顶点，更新度数，以及记录哪些边是不需要插入的。
                int merge_col = col_que[0];
                iter_vtx = merge_col;
                dg->degree[merge_col] = 0;
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                    {
                        int neigh_row = dg->cur_edge[start_ptr];
                        if (dg->degree[neigh_row] && neigh_row != -1) // 更新列顶点的度数
                        {
                            mergearray[neigh_row] = merge_col;
                            dg->degree[merge_col]++;
                        }
                    }

                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (iter_vtx != -1);

                for (int i = 1; i < col_que_end; i++)
                {
                    int delet_col = col_que[i];
                    iter_vtx = delet_col;
                    stack_last = 0;
                    do
                    {
                        for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                        {
                            if (dg->degree[dg->cur_edge[start_ptr]] && dg->cur_edge[start_ptr] != -1) // 邻居顶点存在
                            {
                                if (mergearray[dg->cur_edge[start_ptr]] != col_que[0]) // 是需要记录的行
                                {
                                    mergearray[dg->cur_edge[start_ptr]] = col_que[0];
                                    stack_vtx_neigh[stack_last] = dg->cur_edge[start_ptr];
                                    stack_edge_id[stack_last] = start_ptr;
                                    stack_last++;
                                }
                                else // 已经被记录的行需要被移除，防止连接时记录了两条一样的边
                                {
                                    dg->degree[dg->cur_edge[start_ptr]]--; // 更新被连接的边的度数
                                    if (dg->degree[dg->cur_edge[start_ptr]] == 1)
                                    {
                                        bucket1.push(dg->cur_edge[start_ptr]);
                                    }
                                    else if (dg->degree[dg->cur_edge[start_ptr]] == 2)
                                    {
                                        bucket2.push(dg->cur_edge[start_ptr]);
                                    }

                                    // dg->cur_edge[start_ptr] = -1;
                                    dg->cur_edge[start_ptr] = dg->cur_edge[dg->vtx_ptr_end[iter_vtx] - 1];
                                    dg->origin_edge_end[start_ptr] = dg->origin_edge_end[dg->vtx_ptr_end[iter_vtx] - 1];
                                    dg->origin_edge_start[start_ptr] = dg->origin_edge_start[dg->vtx_ptr_end[iter_vtx] - 1];

                                    dg->vtx_ptr_end[iter_vtx]--;
                                    start_ptr--;
                                }
                            }
                        }

                        iter_vtx = dg->vtx_link_next[iter_vtx]; // 更新迭代顶点
                    } while (iter_vtx != -1);

                    if (stack_last == 0)
                    {
                        delet_vertex(delet_col, dg);
                        continue;
                    }

                    // iter_vtx = col_que[0];
                    // while (dg->vtx_link_next[iter_vtx] != -1)
                    // {
                    //     iter_vtx = dg->vtx_link_next[iter_vtx];
                    // }

                    // if (dg->vtx_ptr_start[iter_vtx + 1] - dg->vtx_ptr_end[iter_vtx] >= stack_last) // 最后一个顶点有问题！！！！！！！
                    // {
                    //     insert_edgelist_to_dynamic_graph(delet_col, stack_vtx_neigh, stack_edge_id, stack_last, col_que[0], dg);
                    // }
                    // else
                    // {
                    connect_edgelist_to_dynamic_graph(delet_col, stack_vtx_neigh, stack_edge_id, stack_last, col_que[0], dg);
                    // }
                }

                // verify_vtx(col_que[0], dg, stack_vtx, verify_array);
                // 按照轮次处理不处理新生成的列顶点
                if (dg->degree[col_que[0]] == 1)
                {
                    bucket1.push(col_que[0]);
                }
                else if (dg->degree[col_que[0]] == 2)
                {
                    bucket2.push(col_que[0]);
                }
            }
        }
        // verify_edge_connection(1212983, 1455456, dg);
    }

    // 记录剩余的顶点
    int num_of_vtx = 0;
    for (int i = 0; i < dg->n; i++)
    {
        if (dg->degree[i] > 0)
        {
            num_of_vtx++;
        }
    }
    cout << "剩下的顶点数量为" << num_of_vtx << endl;

    auto end_kasi = std::chrono::system_clock::now();
    cout << "存储的匹配基数为   " << deg1count + deg2count << " 其中deg1数量为 " << deg1count << " deg2数量为 " << deg2count << endl;

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // of << merge_operations << "    ";
    // of << deg1count + deg2count << "    " << endl;
    // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // of << elapsed_seconds_kasi.count() << endl;
    // of.close();
    // printheapmemory();

    free(mergearray);
    free(row_que);
    free(col_que);
    free(visited[0]);
    free(visited[1]);
    free(visited[2]);
    free(visited);

    free(stack_edge_id);
    free(stack_vtx_neigh);
}

void MVM_raw(graph *og, dynamic_graph *dg, match_tree *mt, int times)
{
    int u, w;
    int start_ptr, end_ptr;
    int start_ptr1, end_ptr1;
    int temp_vtx;
    int temp_vtx_neigh;

    int original_vtx_a, original_vtx_b;

    int iter_vtx;
    stack<int> bucket1;
    stack<int> bucket2;

    // rg->treeedge.resize(og->n);
    // rg->treedge_twin_1.resize(og->n);
    // rg->treedge_twin_2.resize(og->n);
    // rg->treevertexdegree = (int *)calloc(og->n, sizeof(int));

    // 初始化匹配

    int *mergearray = (int *)malloc(sizeof(int) * og->n);
    memset(mergearray, -1, sizeof(int) * og->n);

    // int *verify_array = (int *)malloc(sizeof(int) * og->n);
    // memset(verify_array, -1, sizeof(int) * og->n);

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // auto start_kasi = std::chrono::system_clock::now();
    init_dynamic_graph(og, dg, bucket1, bucket2, times);
    init_match_tree(mt, og);
    // auto end_kasi = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // of << elapsed_seconds_kasi.count() << "    ";
    // of.close();

    // 存储搜索图，对于搜索图中顶点使用bfs，因此使用记录
    int *row_que = (int *)malloc(sizeof(int) * og->nrows);
    unordered_set<int> col_que_mirror;
    int *col_que = (int *)malloc(sizeof(int) * og->nrows);
    int row_que_start, row_que_end, col_que_start, col_que_end;
    // 记录访问顶点               1：记录是否被处理过（记录起始顶点），2：记录被访问次数，3：记录是否已经被加入子图
    int **visited = (int **)malloc(sizeof(int *) * 3);
    visited[0] = (int *)malloc(sizeof(int) * og->n);
    visited[1] = (int *)malloc(sizeof(int) * og->n);
    visited[2] = (int *)malloc(sizeof(int) * og->n);
    memset(visited[0], -1, sizeof(int) * og->n);
    memset(visited[1], -1, sizeof(int) * og->n);
    memset(visited[2], -1, sizeof(int) * og->n);

    // 存储需要移除的列顶点的邻居
    int *stack_vtx_neigh = (int *)malloc(sizeof(int) * og->n);
    int *stack_edge_id = (int *)malloc(sizeof(int) * og->n);
    int stack_last;

    int temp_num = 0;
    int deg1count = 0;
    int deg2count = 0;

    // 两种处理方式：全部顶点：轮次。行顶点：轮次，列顶点：轮次。是否等价。

    auto start_kasi = std::chrono::system_clock::now();

    int merge_operations = 0;
    // 主循环，直到处理完所有的边? 并且度数一和度数二的桶都为空
    while (!bucket1.empty() || !bucket2.empty()) //
    {
        // 度数一的顶点
        while (!bucket1.empty())
        {
            int v = bucket1.top();
            bucket1.pop();
            if (dg->degree[v] == 1) // 隐式的删除方法导致查找邻居需要访问邻接表中所有的元素
            {
                iter_vtx = v;
                u = -1;
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                    {
                        temp_vtx = dg->cur_edge[start_ptr];
                        if (dg->degree[temp_vtx] && temp_vtx != -1) // 是否需要大于0!!!!!!!!!!!!!!!!!!!!
                        {
                            u = temp_vtx;
                            break;
                        }
                    }

                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (u == -1 && iter_vtx != -1); // 找到u直接跳出循环

                // 将该边的原始边加入匹配
                dg->matching[dg->origin_edge_start[start_ptr]] = dg->origin_edge_end[start_ptr];
                dg->matching[dg->origin_edge_end[start_ptr]] = dg->origin_edge_start[start_ptr];

                // 移除顶点v
                delet_vertex(v, dg);

                // 更新顶点u相关的顶点

                iter_vtx = u; // 处理u的所有邻居
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++) // 访问当前顶点的边表
                    {
                        temp_vtx = dg->cur_edge[start_ptr];
                        // 更新tempvtx的状态，如何将tempvtx的邻居边表中u移除？ 只能隐式的移除，否则会增加时间复杂度
                        if (dg->degree[temp_vtx] && temp_vtx != -1) // 度数是准的，顶点指针不一定准
                        {
                            dg->degree[temp_vtx]--;
                            if (dg->degree[temp_vtx] == 1)
                                bucket1.push(temp_vtx);
                            else if (dg->degree[temp_vtx] == 2)
                                bucket2.push(temp_vtx);
                        }
                    }

                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (iter_vtx != -1);

                // 移除顶点u
                delet_vertex(u, dg);

                temp_num++;
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

            if (dg->degree[v] == 2)
            {
                merge_operations++;
                //==================================================================================================================================================//
                //==============================================================搜索==================================================================================//
                //==================================================================================================================================================//
                row_que_start = 0, row_que_end = 0, col_que_start = 0, col_que_end = 0;

                row_que[row_que_end++] = v;
                visited[2][v] = v;

                // 查找顶点v的两个邻居顶点
                int edge_id[2];
                iter_vtx = v;
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                    {
                        temp_vtx = dg->cur_edge[start_ptr];
                        if (dg->degree[temp_vtx] && temp_vtx != -1)
                        {
                            edge_id[col_que_end] = start_ptr;
                            col_que[col_que_end++] = temp_vtx;
                            col_que_mirror.insert(temp_vtx);
                            visited[2][temp_vtx] = v;
                        }
                    }
                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (col_que_end < 2 && iter_vtx != -1);

                augment_recover_graph(edge_id[0], mt, dg);
                augment_recover_graph(edge_id[1], mt, dg);
                // add_twins(e1, e2, rg); // 每次处理一个pair直接构建图

                while (col_que_start != col_que_end)
                {
                    int quecol = col_que[col_que_start++];

                    // 遍历quecol的每一个邻居
                    iter_vtx = quecol;
                    do
                    {
                        for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                        {
                            int temp_row = dg->cur_edge[start_ptr];

                            if (dg->degree[temp_row] == 0 || visited[2][temp_row] == v || temp_row == -1) // 是否要判断度数                             访问的每一个行顶点要判断轮次
                            {
                                continue;
                            }
                            edge_id[0] = start_ptr; // pair的第一条边
                            // 对于可以处理的顶点，更新其起始顶点
                            if (visited[0][temp_row] != v)
                            {
                                visited[0][temp_row] = v;
                            }

                            // 判断行顶点是否满足要求，直接使用差集操作进行判断
                            int sub_num = 0, record_col, record_ptr;
                            // 处理度数存在的列，并且该列要不等于-1，并且没有被加入过子图
                            int iter_vtx1 = temp_row;
                            do
                            {
                                for (start_ptr1 = dg->vtx_ptr_start[iter_vtx1]; start_ptr1 < dg->vtx_ptr_end[iter_vtx1]; start_ptr1++)
                                {
                                    int temp_col = dg->cur_edge[start_ptr1];

                                    // 没有被删除的，没有加入过的
                                    if (dg->degree[temp_col] && temp_col != -1 && visited[2][temp_col] != v)
                                    {
                                        // 如果不在已有集合里面，差集加一，并且对其进行更新
                                        if (col_que_mirror.find(temp_col) != col_que_mirror.end())
                                        {
                                            sub_num++;
                                            record_col = temp_col;
                                            record_ptr = start_ptr1;
                                        }
                                    }
                                }
                                iter_vtx1 = dg->vtx_link_next[iter_vtx1];
                            } while (sub_num <= 1 && iter_vtx1 != -1);

                            if (sub_num == 0 || sub_num == 1)
                            {
                                row_que[row_que_end++] = temp_row;
                                visited[2][temp_row] = v;

                                if (sub_num == 1)
                                {
                                    col_que[col_que_end++] = record_col;
                                    col_que_mirror.insert(record_col);
                                    visited[2][record_col] = v;

                                    // 插入匹配边需要信息
                                    augment_recover_graph(start_ptr, mt, dg);
                                    augment_recover_graph(record_ptr, mt, dg);
                                }
                            }
                            if (row_que_end == col_que_end)
                            {
                                break;
                            }
                        }

                        if (row_que_end == col_que_end)
                        {
                            break;
                        }

                        iter_vtx = dg->vtx_link_next[iter_vtx];
                    } while (iter_vtx != -1);

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
                    delet_vertex(row_que[i], dg);
                }

                // 处理列顶点
                // 将度数最高的列顶点放在首位
                int max_vtx_idx = 0;
                for (int i = 1; i < col_que_end; i++)
                {
                    if (dg->degree[col_que[i]] > dg->degree[col_que[max_vtx_idx]])
                    {
                        max_vtx_idx = i;
                    }
                }
                swap(col_que[0], col_que[max_vtx_idx]);

                // 遍历具有最高度数的列顶点，更新度数，以及记录哪些边是不需要插入的。
                int merge_col = col_que[0];
                iter_vtx = merge_col;
                dg->degree[merge_col] = 0;
                do
                {
                    for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                    {
                        int neigh_row = dg->cur_edge[start_ptr];
                        if (dg->degree[neigh_row] && neigh_row != -1) // 更新列顶点的度数
                        {
                            mergearray[neigh_row] = merge_col;
                            dg->degree[merge_col]++;
                        }
                    }

                    iter_vtx = dg->vtx_link_next[iter_vtx];
                } while (iter_vtx != -1);

                for (int i = 1; i < col_que_end; i++)
                {
                    int delet_col = col_que[i];
                    iter_vtx = delet_col;
                    stack_last = 0;
                    do
                    {
                        for (start_ptr = dg->vtx_ptr_start[iter_vtx]; start_ptr < dg->vtx_ptr_end[iter_vtx]; start_ptr++)
                        {
                            if (dg->degree[dg->cur_edge[start_ptr]] && dg->cur_edge[start_ptr] != -1) // 邻居顶点存在
                            {
                                if (mergearray[dg->cur_edge[start_ptr]] != col_que[0]) // 是需要记录的行
                                {
                                    mergearray[dg->cur_edge[start_ptr]] = col_que[0];
                                    stack_vtx_neigh[stack_last] = dg->cur_edge[start_ptr];
                                    stack_edge_id[stack_last] = start_ptr;
                                    stack_last++;
                                }
                                else // 已经被记录的行需要被移除，防止连接时记录了两条一样的边
                                {
                                    dg->degree[dg->cur_edge[start_ptr]]--; // 更新被连接的边的度数
                                    if (dg->degree[dg->cur_edge[start_ptr]] == 1)
                                    {
                                        bucket1.push(dg->cur_edge[start_ptr]);
                                    }
                                    else if (dg->degree[dg->cur_edge[start_ptr]] == 2)
                                    {
                                        bucket2.push(dg->cur_edge[start_ptr]);
                                    }

                                    // dg->cur_edge[start_ptr] = -1;
                                    dg->cur_edge[start_ptr] = dg->cur_edge[dg->vtx_ptr_end[iter_vtx] - 1];
                                    dg->origin_edge_end[start_ptr] = dg->origin_edge_end[dg->vtx_ptr_end[iter_vtx] - 1];
                                    dg->origin_edge_start[start_ptr] = dg->origin_edge_start[dg->vtx_ptr_end[iter_vtx] - 1];

                                    dg->vtx_ptr_end[iter_vtx]--;
                                    start_ptr--;
                                }
                            }
                        }

                        iter_vtx = dg->vtx_link_next[iter_vtx]; // 更新迭代顶点
                    } while (iter_vtx != -1);

                    if (stack_last == 0)
                    {
                        delet_vertex(delet_col, dg);
                        continue;
                    }

                    connect_edgelist_to_dynamic_graph(delet_col, stack_vtx_neigh, stack_edge_id, stack_last, col_que[0], dg);
                }

                // verify_vtx(col_que[0], dg, stack_vtx, verify_array);
                // 按照轮次处理不处理新生成的列顶点
                if (dg->degree[col_que[0]] == 1)
                {
                    bucket1.push(col_que[0]);
                }
                else if (dg->degree[col_que[0]] == 2)
                {
                    bucket2.push(col_que[0]);
                }

                col_que_mirror.clear();
            }
        }
        // verify_edge_connection(1212983, 1455456, dg);
    }

    // 记录剩余的顶点
    int num_of_vtx = 0;
    for (int i = 0; i < dg->n; i++)
    {
        if (dg->degree[i] > 0)
        {
            num_of_vtx++;
        }
    }
    cout << "剩下的顶点数量为" << num_of_vtx << endl;

    auto end_kasi = std::chrono::system_clock::now();
    cout << "存储的匹配基数为   " << deg1count + deg2count << " 其中deg1数量为 " << deg1count << " deg2数量为 " << deg2count << endl;

    // ofstream of;
    // of.open("TIME.txt", ios ::app);
    // of << merge_operations << "    ";
    // of << deg1count + deg2count << "    " << endl;
    // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // of << elapsed_seconds_kasi.count() << endl;
    // of.close();
    // printheapmemory();

    free(mergearray);
    free(row_que);
    free(col_que);
    free(visited[0]);
    free(visited[1]);
    free(visited[2]);
    free(visited);

    free(stack_edge_id);
    free(stack_vtx_neigh);
}
