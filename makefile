CXX=g++
CC=gcc



CFLAGS_CUSTOM= -O9000  -std=c++11
CFLAGS_LZ4+= -w -Wall -std=gnu99 -DUSE_THREADS  -fstrict-aliasing -Iext $(DEFS)
CFLAGS_BLIGHT+= -DNDEBUG -Ofast -flto -march=native -mtune=native -g -std=c++11 -pipe -lz -fopenmp -msse4 -Ilz4
INC=blight.h bbhash.h common.h
LZ4H=lz4/lz4frame.o lz4/lz4.o lz4/xxhash.o lz4/lz4hc.o



#HERE snippet should be your program using Blight
snippet: snippet.o blight.o utils.o $(LZ4H)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

snippet.o: snippet.cpp
	$(CXX) -o $@ -c $< $(CFLAGS_CUSTOM)



#Here snippet should also be your program
EXEC= snippet bench_blight
#bench_blight is only here for benchmark purpose and can be removed



bench_blight: bench_blight.o blight.o utils.o $(LZ4H)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

bench_blight.o: bench_blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

utils.o: utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

blight.o: blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)
	
lz4/lz4frame.o: lz4/lz4frame.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)

lz4/lz4.o: lz4/lz4.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)

lz4/lz4hc.o: lz4/lz4hc.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)
	
lz4/xxhash.o: lz4/xxhash.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4) 



all: $(EXEC)



clean:
	rm -rf *.o
	rm -rf $(EXEC)
	


rebuild: clean $(EXEC)
