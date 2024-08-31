#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <random>
#include <fstream>
#include <string.h>
#include <sstream>
#include <malloc.h>
#include <cstdio>
#include <string>
#include <stdlib.h>
#include <sys/resource.h>
#include "../include/graphgenBP.h"
#include "../include/kasi.h"
#include "../include/dgmvm.h"

#include <dirent.h> //读取文件相关
#include <iomanip>
#include <numeric>
#include <fstream>

using namespace std;

void getFiles(string path, vector<string> &filenames)
{
    DIR *pDir;
    struct dirent *ptr;
    if (!(pDir = opendir(path.c_str())))
    {
        cout << "Folder doesn't Exist!" << endl;
        return;
    }
    while ((ptr = readdir(pDir)) != 0)
    {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0)
        {
            filenames.push_back(path + "/" + ptr->d_name);
        }
    }
    closedir(pDir);
}
void shuffle_perm(int *rand_perm, int nov, int my_seed)
{
    minstd_rand simple_rand;
    simple_rand.seed(my_seed);
    shuffle(&rand_perm[0], &rand_perm[nov - 1], simple_rand);
}
void shuffle_list(int **rand_edge_list, int nnz, int my_seed)
{
    minstd_rand simple_rand;
    simple_rand.seed(my_seed);
    shuffle(&rand_edge_list[0], &rand_edge_list[nnz - 1], simple_rand);
}

void shuffle_graph(graph *og)
{
    // 将所有的边记录在edge_list中
    int nnz = og->m / 2;
    int **rand_edge_list = (int **)malloc(sizeof(int *) * nnz);
    int edge_num = 0;
    for (int i = 0; i < og->nrows; ++i)
    {
        for (int j = og->vtx_pointer[i]; j < og->vtx_pointer[i + 1]; ++j)
        {
            rand_edge_list[edge_num] = (int *)malloc(2 * sizeof(int));
            rand_edge_list[edge_num][0] = i;
            rand_edge_list[edge_num++][1] = og->endV[j];
        }
    }
    int r_nov = og->nrows;
    int c_nov = og->n - og->nrows;
    int *r_perm = (int *)malloc(sizeof(int) * r_nov);
    int *c_perm = (int *)malloc(sizeof(int) * c_nov);
    int *r_degs = (int *)malloc(sizeof(int) * r_nov);
    int *c_degs = (int *)malloc(sizeof(int) * c_nov);

    // 对于每一条边记录两个值
    shuffle_list(rand_edge_list, nnz, 1); // 打乱边表

    for (int i = 0; i < r_nov; ++i)
        r_perm[i] = i;
    for (int i = 0; i < c_nov; ++i)
        c_perm[i] = i;

    shuffle_perm(r_perm, r_nov, 1); // 打乱顶点
    shuffle_perm(c_perm, c_nov, 1);
    // 对度数进行一个排序

    // 获取顶点度数
    for (int i = 0; i < r_nov; ++i)
        r_degs[r_perm[i]] = og->vtx_pointer[i + 1] - og->vtx_pointer[i];
    for (int i = 0; i < c_nov; ++i)
        c_degs[c_perm[i]] = og->vtx_pointer[i + 1 + r_nov] - og->vtx_pointer[i + r_nov];

    // 重新构建指针
    og->vtx_pointer[0] = 0;
    for (int i = 1; i <= r_nov; ++i)
    {
        og->vtx_pointer[i] = og->vtx_pointer[i - 1] + r_degs[i - 1];
        r_degs[i - 1] = og->vtx_pointer[i] - 1; // 重置度数数组为 指向的最后一条边
    }
    og->vtx_pointer[r_nov] = nnz;
    for (int i = 1; i <= c_nov; ++i)
    {
        og->vtx_pointer[i + r_nov] = og->vtx_pointer[i - 1 + r_nov] + c_degs[i - 1];
        c_degs[i - 1] = og->vtx_pointer[i + r_nov] - 1;
    }

    int u, v, ru, rv;
    // 处理每一条边，将边重新放入矩阵
    for (int i = 0; i < nnz; ++i)
    {
        // 取出一条边
        u = rand_edge_list[i][0];         // 行
        v = rand_edge_list[i][1] - r_nov; // 列

        ru = r_perm[u]; // 行现在所在位置
        rv = c_perm[v]; // 列现在所在位置 从0开始

        // 从后往前放
        og->endV[r_degs[ru]--] = rv + r_nov;
        og->endV[c_degs[rv]--] = ru;
    }

    free(r_perm);
    free(c_perm);
    free(r_degs);
    free(c_degs);
    for (int i = 0; i < nnz; i++)
    {
        free(rand_edge_list[i]);
    }
    free(rand_edge_list);
}

void printheapmemory()
{
    // 同一个程序会测出相同的值
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line))
    {
        if (line.substr(0, 6) == "VmHWM:") // VmSize VmRSS
        {
            std::istringstream iss(line.substr(7));
            size_t memory;
            iss >> memory;
            ofstream of;
            of.open("memory_comsum.txt", ios::app);
            of << line.substr(0, 6) << memory << " KB" << endl;
            of.close();
            std::cout << line.substr(0, 6) << memory << " KB" << std::endl;
        }
        if (line.substr(0, 6) == "VmRSS:") // VmSize VmRSS
        {
            std::istringstream iss(line.substr(7));
            size_t memory;
            iss >> memory;
            ofstream of;
            of.open("memory_comsum.txt", ios::app);
            of << line.substr(0, 6) << memory << " KB" << endl;
            of.close();
            std::cout << line.substr(0, 6) << memory << " KB" << std::endl;
            break;
        }
    }
    // 测出的栈大小不变
    // struct rusage ru;
    // getrusage(RUSAGE_SELF, &ru);
    // printf("Max resident set size (KB): %ld\n", ru.ru_maxrss);
}

void init_graph(graph *og, graph *og1)
{
    og1->n = og->n;
    og1->m = og->m;
    og1->nrows = og->nrows;
    og1->endV = new long[og->m];
    memcpy(og1->endV, og->endV, sizeof(long) * og->m);
    og1->vtx_pointer = new long[og->n + 1];
    memcpy(og1->vtx_pointer, og->vtx_pointer, sizeof(long) * (og->n + 1));
    og1->weight = nullptr;
}

void extract_lines_to_file()
{
    string source_filename = "memory_comsum_tree.txt";
    string destination_filename = "memory_comsum.txt";
    std::ifstream source_file(source_filename);
    std::ofstream destination_file(destination_filename, ios::app);

    if (!source_file.is_open())
    {
        std::cerr << "Error opening source file: " << source_filename << std::endl;
        return;
    }

    if (!destination_file.is_open())
    {
        std::cerr << "Error opening destination file: " << destination_filename << std::endl;
        source_file.close();
        return;
    }

    std::string line;
    int line_count = 0;

    while (std::getline(source_file, line))
    {
        line_count++; // Increment line count
                      // Check if the current line number is in the list of lines to extract
        if (line_count % 5 == 3 || line_count % 5 == 0)
        {
            destination_file << line << std::endl; // Write the line to the destination file
        }
    }

    // Close the files
    source_file.close();
    destination_file.close();
}

int main(int argc, char **argv)
{
    // ceshi();
    extract_lines_to_file();
    if (argc != 3)
    {
        printf("Usage: ./main fileName numThreads\n");
        return -1;
    }

    long numThreads = atoi(argv[2]);
    omp_set_num_threads(numThreads);
#pragma omp parallel
    {
        long nthreads = omp_get_num_threads();
        long nprocs = omp_get_num_procs();
        long tid = omp_get_thread_num();
    }

    vector<string> fil;
    // getFiles("/home/wuguang/bipartite/data/unperfect", fil);
    // for (int iter = 0; iter < 5; iter++)
    // {
    //     for (int filenum = 0; filenum < 30; filenum++)
    //     {
    // printheapmemory();
    // const char *inFile = fil[filenum].c_str();
    const char *inFile = argv[1];
    FILE *fp;
    fp = fopen(inFile, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error! Could not open input file. Exiting ...\n");
        exit(-1);
    }
    fclose(fp);

    graph *og = (graph *)malloc(sizeof(graph));
    graph *og_dg = (graph *)malloc(sizeof(graph));
    graph *og_bg = (graph *)malloc(sizeof(graph));
    graph *og_tg = (graph *)malloc(sizeof(graph));
    graph *og_hg = (graph *)malloc(sizeof(graph));

    fast_mtx_read_build(inFile, og);

    long isolated_rows = 0, isolated_cols = 0;
#pragma omp parallel
    {
        long tisor = 0, tisoc = 0;
#pragma omp for
        for (long u = 0; u < og->nrows; u++)
            if (og->vtx_pointer[u + 1] == og->vtx_pointer[u])
                tisor++;
#pragma omp for
        for (long u = og->nrows; u < og->n; u++)
            if (og->vtx_pointer[u + 1] == og->vtx_pointer[u])
                tisoc++;

        __sync_fetch_and_add(&isolated_rows, tisor); // 原子操作
        __sync_fetch_and_add(&isolated_cols, tisoc);
    }

    printf("===================================\n");
    printf("Problem Statistics\n");
    printf("===================================\n");
    printf("vertices : %ld\n", og->n);
    printf("rows = %ld cols = %ld\n", og->nrows, og->n - og->nrows);
    printf("Number of edges : %ld\n", og->m);
    printf("===================================\n");

    long NV = og->n;
    long nrows = og->nrows;

    shuffle_graph(og);

    // #ifdef TEST

    init_graph(og, og_bg);
    // init_graph(og, og_dg);
    // init_graph(og, og_tg);
    // init_graph(og, og_hg);
    // free_graph(og);
    // ofstream of;
    // of.open("memory_comsum.txt", ios ::app);
    // of << "Begin testing matrix" << inFile << endl;
    // of.close();
    //==================================================================================================================================================//
    //================ Test kernelization algorithms ===================================================================================================//
    //==================================================================================================================================================//

    // directly using the exact algorithm
    // auto start_origin = std::chrono::system_clock::now();
    // long *mateI = (long *)malloc(NV * sizeof(long));
    // memset(mateI, -1, sizeof(long) * NV);
    // Pothen_Fan_Fairnes(og, mateI);
    // free(mateI);
    // auto end_origin = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_origin = end_origin - start_origin;
    // std::cout << "PF time: " << elapsed_seconds_origin.count() << "s\n";

    // // 原始kasi
    // printheapmemory();
    // auto start_kasi = std::chrono::system_clock::now();
    // // of << og->n << "    " << og->m << endl;
    // rec_graph rg;
    // KS2basic(og, &rg);
    // // reduce_the_graph(og, &rg);
    // // of << og->n << "    " << og->m << endl;
    // // long *mate = Pothen_Fan_Fairnes(og, mateI);
    // // recover_the_matching(og, &rg, mate);
    // auto end_kasi = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_kasi = end_kasi - start_kasi;
    // free_rec_graph(&rg, og->n);
    // free_graph(og);
    // std::cout << "KASI time: " << elapsed_seconds_kasi.count() << "s\n";

    // 使用GMVM作为启发时

    // memset(mateI, -1, sizeof(long) * NV);
    // printheapmemory();
    // auto start_dgmvm = std::chrono::system_clock::now();
    // dynamic_graph dg;
    // match_tree mt;
    // DGMVM(og_dg, &dg, &mt, 1);
    // auto end_dgmvm = std::chrono::system_clock::now();
    // // reduce_the_dynamic_graph(og_dg, &dg, &mt);
    // // auto end_dgmvm1 = std::chrono::system_clock::now();
    // // long *mate1 = Pothen_Fan_Fairnes(og_dg, mateI);
    // // auto end_dgmvm2 = std::chrono::system_clock::now();
    // // recover_matching_from_match_tree(og_dg, &dg, &mt, mate1);
    // // auto end_dgmvm3 = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_dgmvm = end_dgmvm - start_dgmvm;
    // // std::chrono::duration<double> elapsed_seconds_dgmvm1 = end_dgmvm1 - end_dgmvm;
    // // std::chrono::duration<double> elapsed_seconds_dgmvm2 = end_dgmvm2 - end_dgmvm1;
    // // std::chrono::duration<double> elapsed_seconds_dgmvm3 = end_dgmvm3 - end_dgmvm2;
    // free_dynamic_graph(&dg); // DGMVM
    // free_match_tree(&mt);
    // free_graph(og_dg);
    // std::cout << "DGMVM time: " << elapsed_seconds_dgmvm.count() << "s\n";
    //==================================================================================================================================================//
    //============================四种间隙的DBMVM算法====================================================================================================//
    //==================================================================================================================================================//
    // memset(mateI, -1, sizeof(long) * NV);theapmemory();
    // printheapmemory();
    auto start_dbmvm1 = std::chrono::system_clock::now();
    dynamic_graph db_g1;
    match_tree db_mt1;
    DBMVM(og_bg, &db_g1, &db_mt1, 1);
    auto end_dbmvm1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_dbmvm1 = end_dbmvm1 - start_dbmvm1;
    free_dynamic_graph(&db_g1); // DBMVM
    free_match_tree(&db_mt1);
    free_graph(og_bg);
    std::cout << "DBMVM time: " << elapsed_seconds_dbmvm1.count() << "s\n";

    // auto start_dbmvm2 = std::chrono::system_clock::now();
    // dynamic_graph db_g2;
    // match_tree db_mt2;
    // DBMVM(og_bg, &db_g2, &db_mt2, 2);
    // auto end_dbmvm2 = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_dbmvm2 = end_dbmvm2 - start_dbmvm2;
    // std::cout << "DBMVM2 time: " << elapsed_seconds_dbmvm2.count() << "s\n";

    // auto start_dbmvm4 = std::chrono::system_clock::now();
    // dynamic_graph db_g4;
    // match_tree db_mt4;
    // DBMVM(og_bg, &db_g4, &db_mt4, 4);
    // auto end_dbmvm4 = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_dbmvm4 = end_dbmvm4 - start_dbmvm4;
    // std::cout << "DBMVM4 time: " << elapsed_seconds_dbmvm4.count() << "s\n";

    // auto start_dbmvm8 = std::chrono::system_clock::now();
    // dynamic_graph db_g8;
    // match_tree db_mt8;
    // DBMVM(og_bg, &db_g8, &db_mt8, 8);
    // auto end_dbmvm8 = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_dbmvm8 = end_dbmvm8 - start_dbmvm8;
    // std::cout << "DBMVM8 time: " << elapsed_seconds_dbmvm8.count() << "s\n";

    // // BMVM
    // auto start_bmvm = std::chrono::system_clock::now();
    // rec_graph bmvm_rg;
    // BMVM(og_bg, &bmvm_rg);
    // auto end_bmvm = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_bmvm = end_bmvm - start_bmvm;
    // std::cout << "BMVM time: " << elapsed_seconds_bmvm.count() << "s\n";

    // of << elapsed_seconds_dbmvm1.count() << "    " << elapsed_seconds_dbmvm2.count() << "    " << elapsed_seconds_dbmvm4.count() << "    " << elapsed_seconds_dbmvm8.count() << "    ";
    //==================================================================================================================================================//
    //=====================================不同的KaSi的效果====================================================================================================//
    //==================================================================================================================================================//
    // // comp_KaSi
    // auto start_comp_kasi = std::chrono::system_clock::now();
    // rec_graph comp_rg;
    // Comp_kaSi(og_bg, &comp_rg);
    // auto end_comp_kasi = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_comp_kasi = end_comp_kasi - start_comp_kasi;
    // std::cout << "Comp_Kasi time: " << elapsed_seconds_comp_kasi.count() << "s\n";

    // auto start_dgmvm = std::chrono::system_clock::now();
    // dynamic_graph dg;
    // match_tree mt;
    // DGMVM(og_dg, &dg, &mt, 1);
    // auto end_dgmvm = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_dgmvm = end_dgmvm - start_dgmvm;
    // hash图

    // printheapmemory();
    // auto start_hash = std::chrono::system_clock::now();
    // hash_graph hg;
    // KaSi_hash(og_hg, &hg);
    // auto end_hash = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_hash = end_hash - start_hash;
    // free_rec_graph(&hg, og_hg->n); // hkasi
    // free_graph(og_hg);
    // std::cout << "hash time: " << elapsed_seconds_hash.count() << "s\n";

    // // // 树图

    // printheapmemory();
    // auto start_tree = std::chrono::system_clock::now();
    // tree_graph tg;
    // KaSi_tree(og_tg, &tg);
    // auto end_tree = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds_tree = end_tree - start_tree;
    // free_rec_graph(&tg, og_tg->n); // TKasi
    // free_graph(og_tg);
    // std::cout << "tree time: " << elapsed_seconds_tree.count() << "s\n";

    // of << elapsed_seconds_kasi.count() << "    " << elapsed_seconds_comp_kasi.count() << "    " << elapsed_seconds_dgmvm.count() << "    " << elapsed_seconds_dbmvm1.count() << "    " << elapsed_seconds_hash.count() << "    " << elapsed_seconds_tree.count() << "    " << endl;
    // of << elapsed_seconds_dgmvm.count() << "    " << elapsed_seconds_dbmvm1.count() << endl;
    // of.close();

    // free_rec_graph(&comp_rg, og->n); // Ckasi
    // free_rec_graph(&bmvm_rg, og->n); // BMVM

    // free_dynamic_graph(&db_g2);
    // free_dynamic_graph(&db_g4);
    // free_dynamic_graph(&db_g8);

    // free_match_tree(&db_mt2);
    // free_match_tree(&db_mt4);
    // free_match_tree(&db_mt8);
    // #endif
    //     }
    // }
    return 0;
}
