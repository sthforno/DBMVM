CC = g++ 
cc = gcc
CFLAGS = -fopenmp -Drestrict=__restrict__ -DNDEBUG -ffast-math -g -pg 
LDFLAGS = -O2

all: main

graphgenBP.o: include/graphgenBP.h source/graphgenBP.cpp include/ThreadedMMReader.h
	$(CC) $(CFLAGS) -c -o graphgenBP.o source/graphgenBP.cpp 

main.o: source/main.cpp include/graphgenBP.h graphgenBP.o
	$(CC) $(CFLAGS) -c -o main.o source/main.cpp 

pf.o: source/pothenfan.cpp include/dgmvm.h include/kasi.h 
	$(CC) $(CFLAGS) -c -o pf.o source/pothenfan.cpp 

mmio.o: source/mmio.c
	$(CC) $(CFLAGS) -Wno-write-strings -c -o mmio.o source/mmio.c

kasi.o: source/kasi_variant.cpp include/kasi.h
	$(CC) $(CFLAGS) -c -o kasi.o source/kasi_variant.cpp

dgmvm.o: source/dgmvm.cpp include/dgmvm.h   
	$(CC) $(CFLAGS) -c -o dgmvm.o source/dgmvm.cpp

dbmvm.o: source/dbmvm.cpp include/dgmvm.h   
	$(CC) $(CFLAGS) -c -o dbmvm.o source/dbmvm.cpp

main: main.o mmio.o kasi.o dgmvm.o dbmvm.o pf.o
	$(CC) $(CFLAGS) $(LDFLAGS)  -o main main.o graphgenBP.o mmio.o kasi.o dgmvm.o dbmvm.o pf.o -lm -lstdc++







clean:
	-rm -f main *.o