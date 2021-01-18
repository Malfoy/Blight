CXX:=g++
CC:=gcc

PREFIX:=/usr/local
LIB:=$(PREFIX)/lib
INC:=$(PREFIX)/include

#Here snippet should be your program
all: libs snippet_static snippet_dynamic bench_blight
#bench_blight is only here for benchmark purpose and can be removed

libs:libBlight.a libBlight.so

#Your compilation flags
CFLAGS_CUSTOM= -O9000  -std=c++11

#Other compilation flags
CFLAGS:=
CFLAGS_BLIGHT+=$(CFLAGS) -DNDEBUG -Ofast -flto -march=native -mtune=native -g -std=c++11 -pipe -lz -fopenmp -msse4 -Ilz4

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
snippet_static: snippet.o libBlight.a
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT) ./libBlight.a -llz4

# example of dynamic library use
snippet_dynamic: snippet.o libBlight.so
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT) -L. -lBlight -llz4

#Benchmark compilation
bench_blight: bench_blight.o 
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT) -L. -lBlight -llz4

bench_blight.o: bench_blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

#Blight object files (BLO)
utils.o: utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)  -fPIC

blight.o: blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)  -fPIC

install:
	mkdir -m 2775 -p $(LIB)
	mkdir -m 2775 -p $(INC)/include
	install -m 0775 *.so $(LIB)
	install -m 0664 *.a $(LIB)
	install -m 0664 *.h $(INC)
	install -m 0664 include/*.h* $(INC)/include

clean:
	rm -f *.o *.a *.so
	rm -f snippet_* bench_blight blight_index.gz
	
rebuild: clean all
