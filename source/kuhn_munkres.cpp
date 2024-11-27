#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <queue>
#include <time.h>
#include <string.h>
#include "../include/kasi.h"
#include "../include/graphgenBP.h"

struct kuhn_munkres
{
    int ncol;
    int nrow;
    int round;

    int record_value;

    int *lx;
    int *ly;

    int *slack;

    int *match; // 列顶点的匹配数组
    int *row_match;

    int *prev; // 存储上一行顶点

    int *weight;

    int *S;
    int *T;

    int bucket_size;
    int *bucket;

    int *tree_row;
    int tree_row_size;
    int *tree_col;
    int tree_col_size;
};

void free_km(kuhn_munkres *km)
{
    free(km->lx);
    free(km->ly);
    free(km->slack);
    free(km->match);
    free(km->row_match);
    free(km->prev);

    free(km->weight);
    free(km->S);
    free(km->T);

    free(km->bucket);   // 存储还没有匹配并且没有处理的行顶点
    free(km->tree_row); // 访问到的在相等子图中的行顶点
    free(km->tree_col); // 访问到的在相等子图中的列顶点
}
void init_km(graph *og, kuhn_munkres *km)
{
    // 构建左右顶标及赋予其初始值
    km->record_value = 0;
    km->nrow = og->nrows;
    km->ncol = og->n - og->nrows;
    if (km->ncol < km->nrow) // 如果列小，直接补齐，从行顶点出发搜索，行小的话总是可以找到增广路径不用改
    {
        km->ncol = km->nrow;
    }

    km->round = 0;

    // 构建左右顶标和松弛值
    km->lx = (int *)malloc(sizeof(int) * km->nrow);
    memset(km->lx, 0, sizeof(int) * km->nrow);
    km->ly = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->ly, 0, sizeof(int) * km->ncol);
    km->slack = (int *)malloc(sizeof(int) * km->ncol);

    // 构建两个匹配数组
    km->match = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->match, -1, sizeof(int) * km->ncol);
    km->row_match = (int *)malloc(sizeof(int) * km->nrow);
    memset(km->row_match, -1, sizeof(int) * km->nrow);

    km->prev = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->prev, -1, sizeof(int) * km->ncol);

    // 构建权重矩阵
    srand((unsigned int)time(NULL));
    km->weight = (int *)malloc(sizeof(int) * km->nrow * km->ncol);
    memset(km->weight, 0, sizeof(int) * km->nrow * km->ncol);
    for (int i = 0; i < km->nrow; i++)
    {
        for (int ptr = og->vtx_pointer[i]; ptr < og->vtx_pointer[i + 1]; ptr++)
        {
            int temp_col = og->endV[ptr] - og->nrows; // 将列顶点的起始位置归为0
            // km->weight[i * km->ncol + temp_col] = int(rand() % 1000 + 1);
            km->weight[i * km->ncol + temp_col] = og->weight[ptr];
        }
    }

    // 初始化行标号lx为各行的最大元素
    for (int i = 0; i < km->nrow; ++i)
    {
        for (int j = 0; j < km->ncol; ++j)
        {
            km->lx[i] = std::max(km->lx[i], km->weight[i * km->ncol + j]);
        }
    }

    // 构建相等子图的访问数组
    km->S = (int *)malloc(sizeof(int) * km->nrow);
    km->T = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->S, -1, sizeof(int) * km->nrow);
    memset(km->T, -1, sizeof(int) * km->ncol);

    // 存储未匹配行顶点
    km->bucket = (int *)malloc(sizeof(int) * km->nrow);
    km->bucket_size = 0;

    // 记录访问到的顶点
    km->tree_col = (int *)malloc(sizeof(int) * km->ncol);
    km->tree_row = (int *)malloc(sizeof(int) * km->nrow);
}

void init_km_with_WU(graph *og, kuhn_munkres *km)
{
    // 先构建km矩阵 （权值矩阵）
    km->nrow = og->nrows;
    km->ncol = og->n - og->nrows;
    if (km->ncol < km->nrow) // 如果列小，直接补齐，从行顶点出发搜索，行小的话总是可以找到增广路径不用改
    {
        km->ncol = km->nrow;
    }
    srand((unsigned int)time(NULL));
    int *weight = (int *)malloc(sizeof(int) * km->nrow * km->ncol);
    memset(weight, 0, sizeof(int) * km->nrow * km->ncol);
    for (int i = 0; i < km->nrow; i++)
    {
        for (int ptr = og->vtx_pointer[i]; ptr < og->vtx_pointer[i + 1]; ptr++)
        {
            int temp_col = og->endV[ptr] - og->nrows; // 将列顶点的起始位置归为0
            weight[i * km->ncol + temp_col] = int(rand() % 1000 + 1);
            // weight[i * km->ncol + temp_col] = og->weight[ptr];
        }
    }
    // 在km矩阵上更新权值

    og->degree = (long *)malloc(sizeof(long) * og->n);
    for (int i = 0; i < og->n; i++)
    {
        og->degree[i] = og->vtx_pointer[i + 1] - og->vtx_pointer[i];
    }

    stack<int> bucket1, bucket2;
    km->record_value = 0;
    for (int i = 0; i < og->nrows; i++)
    {
        if (og->degree[i] == 1)
        {
            bucket1.push(i);
        }
    }
    for (int i = og->nrows; i < og->n; i++)
    {
        if (og->degree[i] == 1)
        {
            bucket2.push(i);
        }
    }
    while (bucket1.size() || bucket2.size())
    {
        while (bucket1.size())
        {
            int temp_row = bucket1.top();
            bucket1.pop();
            if (og->degree[temp_row] == 1)
            {
                og->degree[temp_row] = 0;
                for (int ptr = og->vtx_pointer[temp_row]; ptr < og->vtx_pointer[temp_row + 1]; ptr++)
                {
                    int temp_col = og->endV[ptr];
                    if (og->degree[temp_col] != 0)
                    {
                        int match_value = weight[temp_row * km->ncol + temp_col - og->nrows];
                        if (match_value == 0)
                        {
                            continue;
                        }
                        km->record_value += match_value;
                        weight[temp_row * km->ncol + temp_col - og->nrows] = 0;
                        og->degree[temp_col]--;
                        if (og->degree[temp_col] == 1)
                        {
                            bucket2.push(temp_col);
                        }
                        for (int ptr = og->vtx_pointer[temp_col]; ptr < og->vtx_pointer[temp_col + 1]; ptr++)
                        {
                            int temp_row_2 = og->endV[ptr];
                            if (og->degree[temp_row_2] != 0 && weight[temp_row_2 * km->ncol + temp_col - og->nrows] > 0) // 有值才能够处理
                            {
                                weight[temp_row_2 * km->ncol + temp_col - og->nrows] -= match_value; // 更新相连边的权值

                                if (weight[temp_row_2 * km->ncol + temp_col - og->nrows] <= 0)
                                {
                                    weight[temp_row_2 * km->ncol + temp_col - og->nrows] = 0;
                                    og->degree[temp_col]--;
                                    og->degree[temp_row_2]--; // 是否需要更新度数
                                    if (og->degree[temp_row_2] == 1)
                                    {
                                        bucket1.push(temp_row_2);
                                    }
                                    if (og->degree[temp_col] == 1)
                                    {
                                        bucket2.push(temp_col);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        while (bucket2.size())
        {
            int temp_col = bucket2.top();
            bucket2.pop();
            if (og->degree[temp_col] == 1)
            {
                og->degree[temp_col] = 0;
                for (int ptr = og->vtx_pointer[temp_col]; ptr < og->vtx_pointer[temp_col + 1]; ptr++)
                {
                    int temp_row = og->endV[ptr];
                    if (og->degree[temp_row] != 0)
                    {
                        int match_value = weight[temp_row * km->ncol + temp_col - og->nrows];
                        if (match_value == 0)
                        {
                            continue;
                        }
                        km->record_value += match_value;
                        weight[temp_row * km->ncol + temp_col - og->nrows] = 0;
                        og->degree[temp_row]--;
                        if (og->degree[temp_row] == 1)
                        {
                            bucket1.push(temp_row);
                        }

                        for (int ptr = og->vtx_pointer[temp_row]; ptr < og->vtx_pointer[temp_row + 1]; ptr++)
                        {
                            int temp_col_2 = og->endV[ptr];
                            if (og->degree[temp_col_2] != 0 && weight[temp_row * km->ncol + temp_col_2 - og->nrows] > 0)
                            {
                                weight[temp_row * km->ncol + temp_col_2 - og->nrows] -= match_value;
                                if (weight[temp_row * km->ncol + temp_col_2 - og->nrows] <= 0)
                                {
                                    weight[temp_row * km->ncol + temp_col_2 - og->nrows] = 0;
                                    og->degree[temp_row]--;
                                    og->degree[temp_col_2]--;
                                    if (og->degree[temp_col_2] == 1)
                                    {
                                        bucket2.push(temp_col_2);
                                    }
                                    if (og->degree[temp_row] == 1)
                                    {
                                        bucket1.push(temp_row);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 更新km矩阵
    int old_nrow = km->nrow;
    int old_ncol = km->ncol;

    km->nrow = 0;
    km->ncol = 0;
    int *procject = (int *)malloc(sizeof(int) * og->n);
    memset(procject, -1, sizeof(int) * og->n);
    for (int i = 0; i < og->nrows; i++)
    {
        if (og->degree[i] > 0)
        {
            procject[i] = km->nrow;
            km->nrow++;
        }
    }
    for (int i = og->nrows; i < og->n; i++)
    {
        if (og->degree[i] > 0)
        {
            procject[i] = km->ncol; // 增加列转换为实际列
            km->ncol++;
        }
    }

    if (km->ncol < km->nrow)
    {
        km->ncol = km->nrow;
    }

    km->weight = (int *)malloc(sizeof(int) * km->nrow * km->ncol);
    memset(km->weight, 0, sizeof(int) * km->nrow * km->ncol);
    for (int temp_row = 0; temp_row < og->nrows; temp_row++)
    {
        if (og->degree[temp_row] > 0)
        {
            for (int ptr = og->vtx_pointer[temp_row]; ptr < og->vtx_pointer[temp_row + 1]; ptr++)
            {
                int temp_col = og->endV[ptr];
                if (og->degree[temp_col] > 0) // 对于存在的行和列更新权重矩阵
                {
                    km->weight[procject[temp_row] * km->ncol + procject[temp_col]] = weight[temp_row * old_ncol + temp_col - old_nrow];
                }
            }
        }
    }
    free(weight);
    free(procject);
    free(og->degree);

    // 构建左右顶标和松弛值

    km->prev = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->prev, -1, sizeof(int) * km->ncol);
    km->lx = (int *)malloc(sizeof(int) * km->nrow);
    memset(km->lx, 0, sizeof(int) * km->nrow);
    km->ly = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->ly, 0, sizeof(int) * km->ncol);
    km->slack = (int *)malloc(sizeof(int) * km->ncol);

    // 构建两个匹配数组
    km->match = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->match, -1, sizeof(int) * km->ncol);
    km->row_match = (int *)malloc(sizeof(int) * km->nrow);
    memset(km->row_match, -1, sizeof(int) * km->nrow);
    // 初始化行标号lx为各行的最大元素
    for (int i = 0; i < km->nrow; ++i)
    {
        for (int j = 0; j < km->ncol; ++j)
        {
            km->lx[i] = std::max(km->lx[i], km->weight[i * km->ncol + j]);
        }
    }

    // 构建相等子图的访问数组
    km->S = (int *)malloc(sizeof(int) * km->nrow);
    km->T = (int *)malloc(sizeof(int) * km->ncol);
    memset(km->S, -1, sizeof(int) * km->nrow);
    memset(km->T, -1, sizeof(int) * km->ncol);

    // 存储未匹配行顶点
    km->bucket = (int *)malloc(sizeof(int) * km->nrow);
    km->bucket_size = 0;

    // 记录访问到的顶点
    km->tree_col = (int *)malloc(sizeof(int) * km->ncol);
    km->tree_row = (int *)malloc(sizeof(int) * km->nrow);

    km->round = 0;
}

void update_path(kuhn_munkres *km, int endpoint_col)
{
    int temp_row;
    int temp_col;
    int temp_weight;

    while (endpoint_col != -1)
    {
        temp_row = km->prev[endpoint_col];  // 获取列的前一个顶点
        temp_col = km->row_match[temp_row]; // 最后一个行顶点匹配的是一个-1

        km->match[endpoint_col] = temp_row;
        km->row_match[temp_row] = endpoint_col;

        endpoint_col = temp_col;
    }
}

// 求解过程
bool find_path(int start, kuhn_munkres *km)
{

    // 初始化两个访问数组
    std::queue<int> q;
    q.push(start);

    km->tree_row_size = 0;
    km->tree_col_size = 0;

    km->S[start] = km->round;
    km->tree_row[km->tree_row_size++] = start;
    // 直到找到最大匹配为止？
    while (true)
    {
        // 对于某个行顶点的搜索过程，可以理解，如果不能够访问，记录其能够访问的松弛值，如果可以访问，更新增广路径返回
        while (!q.empty())
        {
            int bound_row = q.front();
            q.pop();

            // 取出一个能够访问的邻居顶点
            for (int bound_col = 0; bound_col < km->ncol; ++bound_col)
            {
                // 如果已经被加入相等子图，则跳过
                if (km->T[bound_col] == km->round)
                    continue;

                // 根据左右顶标判断是否可以访问
                int delta = km->lx[bound_row] + km->ly[bound_col] - km->weight[bound_row * km->ncol + bound_col];

                // 如果存在差值
                if (delta != 0)
                {
                    // 如果松弛值大于差值，更新松弛值
                    if (km->slack[bound_col] > delta)
                    {
                        km->slack[bound_col] = delta; // 可能被重复更新？？？？？？？？
                    }
                }
                else
                {
                    // 将没有加入子图的顶点加入子图，并且设置其前驱
                    km->T[bound_col] = km->round;
                    km->tree_col[km->tree_col_size++] = bound_col;
                    km->prev[bound_col] = bound_row;
                    // 如果列顶点没有匹配？
                    if (km->match[bound_col] == -1)
                    {
                        update_path(km, bound_col);
                        km->round++;
                        return true;
                    }
                    else
                    { // 否则将访问到的顶点的匹配顶点压入路径，将该顶点加入相等子图
                        q.push(km->match[bound_col]);
                        km->S[km->match[bound_col]] = km->round;
                        km->tree_row[km->tree_row_size++] = km->match[bound_col];
                    }
                }
            }
        }

        //==================================================================================================================================================//
        //======================更新顶标,如果没有办法更新顶标，直接返回====================================================================================================================//
        //==================================================================================================================================================//
        int d = 100;
        // 获取没有访问过的顶点的松弛值（松弛值是能够让列顶点的标签降为0的值），构建新的可以访问的顶点
        // for (int j = 0; j < km->ncol; ++j)
        // {
        //     if (km->T[j] != km->round && km->slack[j] < d)
        //         d = km->slack[j];
        // }

        // 只有更新过slack值的列顶点中可能产生最小的d值
        for (int temp_col = 0; temp_col < km->ncol; temp_col++)
        {
            if (km->slack[temp_col] < d)
            {
                d = km->slack[temp_col];
            }
        }
        if (d == 0)
        {
            km->round++;
            return true;
        }

        // 对于访问过的左侧顶点更新顶标
        for (int i = 0; i < km->tree_row_size; i++)
        {
            km->lx[km->tree_row[i]] -= d;
        }
        km->tree_row_size = 0;
        // for (int i = 0; i < km->nrow; ++i)
        // {
        //     if (km->S[i] == km->round)
        //         km->lx[i] -= d;
        // }
        // 对于访问过的右侧顶点更新顶标，对于没有访问过的顶点更新其松弛值（下次只需要更少的值即可构建可以访问的顶点）
        for (int i = 0; i < km->tree_col_size; i++)
        {
            km->ly[km->tree_col[i]] += d;
        }
        km->tree_col_size = 0;

        //==================================================================================================================================================//
        //==========更新访问数组以及初始化路径========================================================================================================================//
        //==================================================================================================================================================//
        km->round++;
        q.push(start);
        km->S[start] = km->round;
        km->tree_row[km->tree_row_size++] = start;
    }
}

// 初始化匹配
void init_match(kuhn_munkres *km)
{
    for (int temp_row = 0; temp_row < km->nrow; ++temp_row)
    {
        for (int temp_col = 0; temp_col < km->ncol; ++temp_col)
        {
            if (km->match[temp_col] == -1 && km->lx[temp_row] + km->ly[temp_col] == km->weight[temp_row * km->ncol + temp_col])
            {
                km->match[temp_col] = temp_row;
                km->row_match[temp_row] = temp_col;
                break;
            }
        }
        if (km->row_match[temp_row] == -1)
        {
            km->bucket[km->bucket_size++] = temp_row;
        }
    }
}

void no_init_match(kuhn_munkres *km)
{
    for (int i = 0; i < km->nrow; ++i)
    {
        km->bucket[km->bucket_size++] = i;
    }
}

void solve(kuhn_munkres *km)
{
    init_match(km);
    // no_init_match(km);
    for (int i = 0; i < km->bucket_size; ++i)
    {
        int process_row = km->bucket[i];
        // 求解每一个顶点都更新一次松弛值
        for (int j = 0; j < km->ncol; ++j)
        {
            km->slack[j] = 100;
        }

        find_path(process_row, km);
    }
}

// 寻找增广路径

int test_kuhn_munkres(graph *og)
{
    // 测试数据
    printheapmemory();
    kuhn_munkres km;

    // init_km(og, &km);
    init_km_with_WU(og, &km);

    // solve(&km);

    // std::cout << "获取到的最大匹配值为:" << std::endl;
    // int find_match_value = 0;
    // for (int temp_row = 0; temp_row < km.nrow; ++temp_row)
    // {
    //     find_match_value += km.weight[temp_row * km.ncol + km.row_match[temp_row]];
    //     cout << "行 " << temp_row << "列 " << km.row_match[temp_row] << "匹配值 " << km.weight[temp_row * km.ncol + km.row_match[temp_row]] << endl;
    // }
    // cout << find_match_value + km.record_value << endl;

    // printheapmemory();
    free_km(&km);
}