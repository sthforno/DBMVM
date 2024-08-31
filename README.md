## DBMVM
kernelization algorithm for bipartite graph matching.

The algorithms we designed are included in dgmvm.cpp, dbmvm.cpp and dgmvm.h. The variants of the kasi algorithm are included in kasi_variant.cpp. OpenMP is used for fast file reading.

### Test
```
make
./main com-LiveJournal.mtx 1
```

### Download the real-life dateset
https://sparse.tamu.edu/

### Construct the synthetic graph in the following format:
```
%%MatrixMarket matrix coordinate integer general

6 6 8
1 1 1
2 2 1
2 3 1
3 4 1
4 5 1
5 3 1
6 5 1
6 2 1
```

The third line contains three numbers representing the number of vertices in each of the two sets of a bipartite graph, followed by the total number of edges. The subsequent lines provide data on the edges, with each line listing the starting vertex and ending vertex IDs