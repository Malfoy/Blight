# Blight

## de Bruijn graph-based index with light memory usage

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


## Bird-eye view:
Blight is an associative data structure, similar to a hashtable, able to index a kmer set.
It associates each kmer with a unique identifier and can identify alien kmers (kmer absent from the index) associated with the identifier -1.
The kmer identifiers are in\[0, N\[ where N is the number of kmer in the index.
Hence, one can see Blight as a Minimal Perfect Hash Function (MPHF) that handle alien kmers.


## Key properties:
A Blight index is
- Deterministic, it produces no false positive or false negative.
- Built from a compacted de Bruijn graph
- Static, once constructed, the indexed kmer set cannot be modified
- Fast and memory-efficient even for the most extensive kmer set


## Graph construction
To be constructed, a Blight index needs a Fasta file, whose sequence contains the kmer to index with no duplicate.
An efficient way to build such a file is to construct a compacted de Bruijn graph from interest sequences.
We recommend using BCALM2 (github.com/GATB/bcalm) to construct such a graph from a set of sequences in Fasta, Fastq, gzipped or not.

The bcalm2 command should look like this:

```sh
bcalm -in my_fasta_file.fa -kmer-size 31 -abundance-min 1 -out my_graph
```

Note that by default, bcalm2 will remove unique kmer if the parameter abundance-min is not set.

## Compilation
After the simple installation using:

```sh
git clone --depth 1  https://github.com/Malfoy/Blight
make -j 4
```

We get a snippet binary that plays with the Blight API and a benchmark binary bench_blight.

We can test the blight index's performance on our system using bench_blight on a compacted de Bruijn graph.

We can start to use the API by modifying the snippet.cpp file.

The makefile shows how to compile Blight and lz4 object files and link them to a program.


## Main parameters
We can tune the Blight index with several parameters impacting its overall performances.

### Kmer_size k
All the words of length k present in the fasta file will be indexed.

### Core_Number c
This parameter defines the number of threads used by Blight to perform queries and construct its index.
The default value is one thread.

### Minimizer size m
The minimizer of a kmer is its minimal word of length m  (according to a hash function)
Blight will use 4^m partitions to index the kmer according to their minimizer.
A larger minimizer size reduces the memory fingerprint up to the point where the constant overhead of O(4^m) is too large.
A larger minimizer size can also accelerate the query time by improving cache coherency.
The default value is 10 for roughly 1M partitions.

### File_number_exponent s
Blight will use 4^s temporary files to construct the index.
Allowing more files can accelerate the construction and diminish the memory used during this step.
The default value is 4 for roughly 256 files.

On most Linux system
```sh
ulimit -n 2000
```

Will allow the use of 2000 open files.
 
### Bits_to_save b
Blight will subsample the kmer position in the index.
This option will save b bits per kmer, but each query will perform at most 2^b verification to find a kmer_Set_Light.
With low b values, the query time's impact is low because of cache mechanisms, but for high b values (b>5), an increment double the expected query time.
The default value is 0 (no subsampling)



## Main functions

The snippet toy program includes all the useful functions. We will describe them in detail here.

### Index construction

A blight index is an object with the given constructor setting the index parameter described previously:
```cpp
kmer_Set_Light(uint64_t kmer_size ,uint64_t core_Number=1, uint64_t minimizer_size=10, uint64_t file_number_exponent=4,  uint64_t bits_to_save=0);
```

Note that another constructor allows to load an index directly from an index file containing a previous index dumped on disk (See index on disk)

```cpp
kmer_Set_Light(const string& index_file);
```


Once created, we can construct the index from a fasta file with the construct_index function.
Blight creates its temporary files in a working directory that MUST exist.

```cpp
	void construct_index(const string& input_file, const string& wdir = "");
```


## Queries

Once constructed, a blight index can perform two kinds of queries from a sequence.
Both functions output a vector with a value for each kmer of the query sequence.
The first element of the vector is the value associated with the leftmost kmer.

Membership queries are perfomed with get_presence_query
```cpp
vector<bool> get_presence_query(const string& seq);
```
True meaning that the kmer is present in the index, False otherwise.

Associative queries with get_hashes_query
```cpp
	vector<int64_t> get_hashes_query(const string& seq);
```
The associated indices are in the range \[0, get_kmer_number()\] for the kmer present in the index, -1 otherwise



## Index on disk

Once constructed, the blight index can be saved on the disk to be reused later with the function dump_disk 

```cpp
void dump_disk(const string& output_file);
```

It can be loaded later with a load constructor
```cpp
kmer_Set_Light(const string& output_file);
```

## Blight as a dictionary
Blight combined with an array of values is similar to a static hashtable.
As show in the snippet, an array of value can be allocated using the function get_kmer_number() that indicate the amount of kmer present in the index.
```cpp
int abundance[blight_index_6.get_kmer_number()]={0};
```
This way, Blight will associate each indexed kmer to a position in the abundance array.


## Iterator
It can be useful to be able to iterate on all the kmer of the index.
This can be done using an iterator, as shown in the snippet.

```cpp
kmer_Set_Light_iterator it(&blight_index);
	
do{
	//We can obtain a binary representation of the kmer as a integer
	kmer kmer_binary(it.get_kmer());
	
	//Or as a sequence
	string kmer_sequence(it.get_kmer_str());
	
	
	//If iterator.next() return false it mean that there is no more kmer in the index
}while(it.next());

```











