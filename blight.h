#ifndef KSL
#define KSL



#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "include/common.h"
#include "include/bbhash.h"
#include "include/bm64.h"
#include "include/encoding.h"
#include "include/zstr.hpp"
#include "include/bmserial.h"
#include "include/robin_hood.h"



#define BMSSE2OPT
#define BMSSE42OPT
#define BMAVX2OPT



#define minimizer_type uint32_t



// Use 64 bits int to represent kmer up to k=32
#define kmer uint64_t



// Use 128 bits int to represent kmer up to k=64
// #define kmer __uint128_t



using namespace std;



typedef boomphf::SingleHashFunctor<kmer> hasher_t;
typedef boomphf::mphf<kmer, hasher_t> MPHF;



struct KmerHasher {
	size_t operator()(const kmer& k) const { return ((uint64_t)k); }
};



struct kmer_context {
	bool isdump;
	vector<pair<uint16_t, uint16_t>> count;
};



struct bucket_minimizer {
	uint32_t skmer_number;
};



struct info_mphf {
	uint64_t mphf_size;
	uint64_t start;
	MPHF* kmer_MPHF;
	uint8_t bit_to_encode;
	bool empty;
};



struct minitig {
	int32_t color;
	uint16_t coverage;
	string sequence;
};



// Represents the cardinality of a pow2 sized set. Allows div/mod arithmetic operations on indexes.
template<typename T>
struct Pow2 {
	Pow2(uint_fast8_t bits)
	  : _bits(bits) {
		assume(bits < CHAR_BIT * sizeof(T), "Pow2(%u > %u)", unsigned(bits), unsigned(CHAR_BIT * sizeof(T)));
	}

	uint_fast8_t bits() const { return _bits; }
	T value() const { return T(1) << _bits; }
	explicit operator T() const { return value(); }
	T max() const { return value() - T(1); }
	Pow2& operator=(const Pow2&) = default;

	Pow2()
	  : _bits(0) {}
	Pow2(const Pow2&) = default;

	friend T operator*(const T& x, const Pow2& y) { return x << y._bits; }
	friend T& operator*=(T& x, const Pow2& y) { return x <<= y._bits; }
	friend T operator/(const T& x, const Pow2& y) { return x >> y._bits; }
	friend T& operator/=(T& x, const Pow2& y) { return x >>= y._bits; }
	friend T operator%(const T& x, const Pow2& y) { return x & y.max(); }
	friend T& operator%=(T& x, const Pow2& y) { return x &= y.max(); }
	Pow2& operator>>=(uint_fast8_t d) {
		_bits -= d;
		return *this;
	}
	Pow2& operator<<=(uint_fast8_t d) {
		_bits += d;
		return *this;
	}
	friend bool operator<(const T& x, const Pow2& y) { return x < y.value(); }
	friend bool operator<=(const T& x, const Pow2& y) { return x < y.value(); }
	friend T operator+(const T& x, const Pow2& y) { return x + y.value(); }
	friend T& operator+=(T& x, const Pow2& y) { return x += y.value(); }
	friend T operator-(const T& x, const Pow2& y) { return x - y.value(); }
	friend T& operator-=(T& x, const Pow2& y) { return x -= y.value(); }

  private:
	uint_fast8_t _bits;
};



class kmer_Set_Light {
  public:
	uint64_t k, m1, m2, m3, minimizer_size_graph;
	uint64_t coreNumber;
	uint64_t bit_saved_sub;
	uint64_t total_nuc_number = 0;
	int color_mode;
	string working_dir;
	double max_divergence_count = 0;

	Pow2<kmer> offsetUpdateAnchor;
	Pow2<kmer> offsetUpdateMinimizer;
	Pow2<uint> mphf_number;
	Pow2<uint> number_superbuckets;
	Pow2<uint> minimizer_number;
	Pow2<uint> minimizer_number_graph;
	Pow2<uint> number_bucket_per_mphf;
	Pow2<uint> bucket_per_superBuckets;
	Pow2<uint> positions_to_check;

	vector<bool> bucketSeq;
	vector<bool> positions;
	bm::bvector<> position_super_kmers = bm::bvector<>(bm::BM_BIT);
	bm::bvector<>::rs_index_type* position_super_kmers_RS;
	mutex positions_mutex[4096];

	bucket_minimizer* all_buckets;
	uint32_t* nuc_minimizer;
	uint64_t* current_pos;
	uint64_t* start_bucket;
	uint32_t* abundance_minimizer_temp;
	uint8_t* abundance_minimizer;
	info_mphf* all_mphf;

	uint64_t number_kmer = 0;
	uint64_t number_super_kmer = 0;
	uint64_t largest_MPHF = 0;
	uint64_t positions_total_size = 0;
	uint64_t positions_int = 0;
	uint64_t read_kmer = 0;
	atomic<uint64_t> number_query;

	uint64_t total_nb_minitigs = 0;
	//~ uint64_t hell_bucket;
	double bit_per_kmer = 0;
	uint32_t largest_bucket_nuc_all = 0;
	uint64_t gammaFactor = 2;
	kmer_Set_Light(uint64_t k_val,uint64_t coreNumber_val=1, uint64_t m1_val=10, uint64_t m3_val=4,  uint64_t bit_to_save=0)
	  : k(k_val)
	  , m1(m1_val)
	  , m2(m1_val)
	  , m3(m3_val)
	  , minimizer_size_graph(m1 + 3)
	  , coreNumber(coreNumber_val)
	  , bit_saved_sub(bit_to_save)

	  , offsetUpdateAnchor(2 * k)
	  , offsetUpdateMinimizer(2 * minimizer_size_graph)
	  , mphf_number(2 * m2)
	  , number_superbuckets(2 * m3)
	  , minimizer_number(2 * m1)
	  , minimizer_number_graph(2 * minimizer_size_graph)
	  , number_bucket_per_mphf(2 * (m1 - m2))
	  , bucket_per_superBuckets(2 * (m1 - m3))
	  , positions_to_check(bit_to_save)
      , position_super_kmers_RS(nullptr) {

		all_buckets = new bucket_minimizer[minimizer_number.value()]();
		all_mphf = new info_mphf[mphf_number.value()];
		for (uint64_t i(0); i < mphf_number; ++i) {
			all_mphf[i].mphf_size = 0;
			all_mphf[i].bit_to_encode = 0;
			all_mphf[i].start = 0;
			all_mphf[i].empty = true;
		}
		for (uint64_t i(0); i < minimizer_number.value(); ++i) {
			all_buckets[i].skmer_number = 0;
		}
		number_query = 0;
	}

	kmer_Set_Light(const string& index_file);

	~kmer_Set_Light() {
        if (position_super_kmers_RS)
    		delete position_super_kmers_RS;
		delete[] all_buckets;

		for (uint64_t i(0); i < mphf_number.value(); ++i) {
			if (not all_mphf[i].empty) {
				delete all_mphf[i].kmer_MPHF;
			}
		}
		delete[] all_mphf;
		struct stat buffer;
		string filerm;
		for (uint64_t i(0); i < number_superbuckets.value(); ++i) {
			filerm = ("_out" + to_string(i));
			if ((stat(filerm.c_str(), &buffer) == 0)) {
				remove(filerm.c_str());
			}
		}
	}

	void read_super_buckets(const string& input_file);
	uint64_t get_kmer_number();
	void create_mphf_mem(uint64_t beg, uint64_t end);
	void create_mphf_disk(uint64_t beg, uint64_t end, bm::bvector<>& position_super_kmers_local);
	void updateK(kmer& min, char nuc);
	void updateRCK(kmer& min, char nuc);
	void updateM(kmer& min, char nuc);
	void updateRCM(kmer& min, char nuc);
	void fill_positions(uint64_t beg, uint64_t end, bm::bvector<>& position_super_kmers_local);
	kmer update_kmer(uint64_t pos, kmer mini, kmer input);
	kmer get_kmer(uint64_t pos, uint64_t mini);
	void print_kmer(kmer num, uint64_t n = 100);
	uint64_t bool_to_int(uint64_t n_bits_to_encode, uint64_t pos, uint64_t start);
	void int_to_bool(uint64_t n_bits_to_encode, uint64_t X, uint64_t pos, uint64_t start);
	kmer update_kmer_local(uint64_t pos, const vector<bool>& V, kmer input);
	string kmer2str(kmer num);
	kmer regular_minimizer(kmer seq);
	void create_super_buckets(const string&);
	void construct_index(const string& input_file, const string& wdir = "");
	vector<kmer> kmer_to_superkmer(const kmer canon, kmer minimizer, int64_t& rank, int64_t& hash);
	int64_t hash_to_rank(const int64_t hash, kmer minimizer);
	int64_t kmer_to_hash(const kmer canon, kmer minimizer);
	string compaction(const string& seq1, const string& seq2, bool);
	void reset();
	void dump_disk(const string& output_file);
	vector<bool> get_presence_query(const string& seq);
	vector<int64_t> get_hashes_query(const string& seq);
	vector<int64_t> get_rank_query(const string& seq);
	void file_query_presence(const string& query_file);
	void file_query_hases(const string& query_file, bool check = true);
	void file_query_rank(const string& query_file);
	void file_query_all_test(const string& query_file, bool);
	__uint128_t rcb(const __uint128_t&);
	uint64_t rcb(const uint64_t&);
	uint64_t canonize(uint64_t x, uint64_t n);
	kmer get_kmer(uint64_t pos);
	void str2bool(const string& str, uint64_t mini);
	void dump_and_destroy(const string& output_file);
	kmer regular_minimizer_pos(kmer seq, uint64_t& position);
	void initialize_buckets();
	//~ //Reindeer
	//~ void write_buffer_count(vector<string>& buffers, zstr::ofstream* out, vector<uint16_t>& headerV, string& seq2dump, int32_t minimi);
	//~ void write_buffer_color(vector<string>& buffers, zstr::ofstream* out, vector<uint8_t>& headerV, string& seq2dump, int32_t minimi);
	uint64_t get_minimizer_from_header(zstr::ifstream& in);
	//~ void merge_super_buckets_mem(const string& input_file, uint64_t number_color, string& out_name,uint64_t number_pass=1, int colormode=1 );
	//~ void get_monocolor_minitigs_mem(vector<robin_hood::unordered_node_map<kmer, kmer_context>>& min2kmer2context,zstr::ofstream* out,const vector<int32_t>& mini,uint64_t number_color, int colormode);
	void read_super_buckets_reindeer(const string& input_file);
};



class kmer_Set_Light_iterator {
public:
	uint64_t rank;
	uint64_t kmer_id;
	uint64_t position;
	uint64_t next_position;
	bm::id64_t rank_position;
	kmer_Set_Light* index_ptr;

	kmer_Set_Light_iterator(kmer_Set_Light* ptdr){
		index_ptr=ptdr;
		kmer_id=rank=1;
		position=0;
		rank_position=0;
		rank_position=(index_ptr->position_super_kmers.get_next(rank_position));
		if (rank_position == 0) {
			next_position = index_ptr->bucketSeq.size() - index_ptr->k;
		}else{
			next_position = rank_position+(index_ptr->k - 1);
		}
	}

	string get_kmer_str()const{
		kmer seqR = index_ptr->get_kmer(position);
		return index_ptr->kmer2str(seqR);
	}


	kmer get_kmer()const{
		return index_ptr->get_kmer(position);
	}

	bool next(){
		kmer_id++;

		if(kmer_id>index_ptr->number_kmer){
			return false;
		}
		position++;
		if(position+index_ptr->k-1>=next_position){
			rank++;
			position=next_position;
			rank_position=(index_ptr->position_super_kmers.get_next(rank_position));
			if (next_position == 0) {
				next_position = index_ptr->bucketSeq.size() - index_ptr->k;
			} else {
				next_position = rank_position + rank  * (index_ptr->k - 1);
			}
		}

		return true;
	}

};



#endif
