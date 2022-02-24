#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <tmmintrin.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include "SubSampler.h"
#include "utils.h"

#include "include/zstr.hpp"




using namespace std;
using namespace chrono;



void Subsampler::updateK(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchor;
}



void Subsampler::updateM(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateMinimizer;
}



void Subsampler::updateRCK(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * k - 2));
}



void Subsampler::updateRCM(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * minimizer_size - 2));
}


uint64_t Subsampler::canonize(uint64_t x, uint64_t n) {
	return min(x, rcbc(x, n));
}


uint64_t Subsampler::regular_minimizer_pos(uint64_t seq, uint64_t& position) {
	uint64_t mini, mmer;
	mmer = seq % minimizer_number;
	mini = mmer        = canonize(mmer, minimizer_size);
	uint64_t hash_mini = (unrevhash(mmer));
	position           = 0;
	for (uint64_t i(1); i <= k - minimizer_size; i++) {
		seq >>= 2;
		mmer          = seq % minimizer_number;
		mmer          = canonize(mmer, minimizer_size);
		uint64_t hash = (unrevhash(mmer));
		if (hash_mini > hash) {
			position  = k - minimizer_size - i;
			mini      = mmer;
			hash_mini = hash;
		}
	}
    return revhash((uint64_t)mini) % minimizer_number;
}


void Subsampler::parse_fasta(const string& input_file) {
	uint64_t total_nuc_number(0);
    uint64_t read_kmer(0);
	istream* input_stream;
	input_stream = new zstr::ifstream(input_file);
	if (not input_stream->good()) {
		cout << "Problem with file opening" << endl;
		exit(1);
	}
	auto out_file=ofstream(("sub"+input_file).c_str(), ofstream::trunc);
	#pragma omp parallel num_threads(coreNumber)
	{
		string ref, useless;
		uint32_t old_minimizer, minimizer;
		while (not input_stream->eof()) {
			ref = useless = "";
			#pragma omp critical(dataupdate)
			{
				getline(*input_stream, useless);
				getline(*input_stream, ref);
				if (ref.size() < k) {
					ref = "";
				} else {
					read_kmer += ref.size() - k + 1;
				}
			}
			// FOREACH UNITIG
			if (not ref.empty() and not useless.empty()) {
				old_minimizer = minimizer = minimizer_number;
				uint64_t last_position(0);
				// FOREACH KMER
				uint64_t seq(str2num(ref.substr(0, k)));
				uint64_t position_min;
				uint64_t min_seq = (str2num(ref.substr(k - minimizer_size, minimizer_size))),
                min_rcseq(rcbc(min_seq, minimizer_size)),
				min_canon(min(min_seq, min_rcseq));
				minimizer         = regular_minimizer_pos(seq, position_min);
				old_minimizer     = minimizer;
				uint64_t hash_min = unrevhash(minimizer);
				uint64_t i(0);
				for (; i + k < ref.size(); ++i) {
					updateK(seq, ref[i + k]);
					updateM(min_seq, ref[i + k]);
					updateRCM(min_rcseq, ref[i + k]);
					min_canon      = (min(min_seq, min_rcseq));
					uint64_t new_h = unrevhash(min_canon);
					// THE NEW mmer is a MINIMIZor
					if (new_h < hash_min) {
						minimizer    = (min_canon);
						hash_min     = new_h;
						position_min = i + k - minimizer_size + 1;
					} else {
						// the previous minimizer is outdated
						if (i >= position_min) {
							minimizer = regular_minimizer_pos(seq, position_min);
							hash_min  = unrevhash(minimizer);
							position_min += (i + 1);
						} else {
						}
					}
					// COMPUTE KMER MINIMIZER
					if (revhash(old_minimizer) % minimizer_number != revhash(minimizer) % minimizer_number) {
						old_minimizer = (revhash(old_minimizer) % minimizer_number);
                        if((i - last_position + 1)==max_superkmer_size){
                            if(old_minimizer%subsampling_rate==0){
                                out_file<<">" + to_string(old_minimizer) + "\n" + ref.substr(last_position, i - last_position + k) + "\n";
                                selected_kmer_number+=(i - last_position + 1);
                                selected_superkmer_number++;
                            }
                        }
						total_kmer_number += (i - last_position + 1);
                        total_superkmer_number++;
						last_position = i + 1;
						old_minimizer = minimizer;
					}
				}
				if (ref.size() - last_position > k - 1) {
					old_minimizer = (revhash(old_minimizer) % minimizer_number);
                    if((ref.size() - last_position + 1)==max_superkmer_size){
                        if(old_minimizer%subsampling_rate==0){
                            out_file<<">" + to_string(old_minimizer) + "\n" + ref.substr(last_position) + "\n";
                            selected_kmer_number+=(ref.size() - last_position + 1);
                            selected_superkmer_number++;
                        }
                    }
                    total_kmer_number += (ref.size() - last_position + 1);
                    total_superkmer_number++;
				}
			}
		}
	}
}



int main(int argc, char** argv) {
	char ch;
	string input, inputfof, query;
	uint k(31);
	uint m1(10);
	uint c(1);
    uint64_t s(8); 

	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:i:")) != -1) {
		switch (ch) {
			case 'i': input = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'm': m1 = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
            case 's': s = stoi(optarg); break;
		}
	}


	if ((input == "" )) {
		cout << "Core arguments:" << endl
		     << "	-i input file" << endl
		     << "	-k kmer size used  (31) " << endl
             << "	-s subsampling used  (8) " << endl
             << "	-m minimize size used  (10) " << endl;
		return 0;
	}else{
        cout<<" I use k="<<k<<" m="<<m1<<" s="<<s<<endl;
        cout<<"Maximal super kmer are of length "<<2*k-m1<<" or "<<k-m1+1<<" kmers" <<endl;
        Subsampler ss(k,m1,s,c);
        ss.parse_fasta(input);
        cout<<"I have seen "<<intToString(ss.total_kmer_number)<<" kmers and I selected "<<intToString(ss.selected_kmer_number)<<" kmers"<<endl;
        cout<<"This mean a pratical subsampling rate of "<<(double)ss.total_kmer_number/ss.selected_kmer_number<<endl;
        cout<<"I have seen "<<intToString(ss.total_superkmer_number)<<" superkmers and I selected "<<intToString(ss.selected_kmer_number)<<" superkmers"<<endl;
        cout<<"This mean a pratical subsampling rate of "<<(double)ss.total_superkmer_number/ss.selected_superkmer_number<<endl;
        cout<<"This mean a mean superkmer size of "<<(double)ss.total_kmer_number/ss.total_superkmer_number<<" kmer per superkmer in the input"<<endl;
        cout<<"This mean a mean superkmer size of "<<(double)ss.selected_kmer_number/ss.selected_superkmer_number<<" kmer per superkmer in the output"<<endl;
    }
	    
}