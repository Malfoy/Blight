CXX=g++
CC=gcc



#Here snippet should be your program
all: libs snippet_static snippet_dynamic bench_blight
#bench_blight is only here for benchmark purpose and can be removed

libs:libBlight.a libBlight.so

#Your compilation flags
CFLAGS_CUSTOM= -O9000  -std=c++11

#Other compilation flags
CFLAGS_BLIGHT+= -DNDEBUG -Ofast -flto -march=native -mtune=native -g -std=c++11 -pipe -lz -fopenmp -msse4 -Ilz4 -fPIC

#Needed object files
BLO=blight.o utils.o

# static library
libBlight.a: $(BLO)
	gcc-ar rvs -o $@ $(BLO)

# dynamic library
libBlight.so: $(BLO)
	$(CC) -shared -o $@ $(BLO)


#Here  you should compile be your program using Blight instead of snippet
snippet: snippet.o $(BLO) 
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT) 

# example of static library use
snippet_static: snippet.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT) libBlight.a -llz4

# example of dynamic library use
snippet_dynamic: snippet.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT) -L. -lBlight -llz4

#Benchmark compilation
bench_blight: bench_blight.o blight.o utils.o 
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT) -L. -lBlight -llz4

bench_blight.o: bench_blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

#Blight object files (BLO)
utils.o: utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

blight.o: blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)
	
clean:
	rm -f *.o *.a *.so
	rm -f snippet_* bench_blight blight_index.gz
	
rebuild: clean all
