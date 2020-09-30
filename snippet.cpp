#include "blight.h"



using namespace std;



int main(int argc, char** argv) {
	


	//SOME VARIABLES TO PLAY WITH
	int kmer_size(31);
	int core_number(31);
	int minimizer_size(10);
	int file_number_exponent(4);
	int subsampling_bits(0);
	string input_fasta_file("My_debruijn_graph.fa");
	string my_query("CTGATCGATCGTACGTAGCTGCTGATCGATCGTACGTACGTACGTCAGT");
	string my_query_bad("ATATATATATATATATATATATATATATATATAT");
	
	
	
	//INDEX INITIZALIZATION
	
	
	
	//Index Initialization with a given kmer size
	kmer_Set_Light blight_index_1(kmer_size);
	
	//Index Initialization allowing the use  of multiple thread for faster construction and queries (default is 1)
	kmer_Set_Light blight_index_2(kmer_size, core_number);
	
	//Index Initialization with a given minimizer size (default is 10)
	kmer_Set_Light blight_index_3(kmer_size, core_number, minimizer_size);
	
	//Index Initialization allowing a custom  amount of temporary file (default is 4 for 256 files)
	kmer_Set_Light blight_index_4(kmer_size, core_number, minimizer_size, file_number_exponent);
	//Blight will create 4^file_number_exponent temporary files (1->4 2>16 3->64 4->256 5->1024)
	
	//Index Initialization using position subsampling to reduce the memory usage of the index (default is 0)
	kmer_Set_Light blight_index_5(kmer_size, core_number, minimizer_size, file_number_exponent, subsampling_bits);
	//A value of N will save up to N bits per kmer but each query can lead up to 2^N comparisons to find a kmer in the index
	
	
	
	//INDEX CONSTRUCTION
	
	
	
	//Construct the index from an input FASTA file, temporary files are put in folder wdir that MUST exist (!)
	blight_index_5.construct_index(input_fasta_file, "wdir");
	//Blight handles .gz and .lz4 files but do not handle FASTQ or multiline FASTA
	
	
	
	//QUERIES
	
	
	
	//Check the presence of all the kmer of the query sequence
	vector<bool> presence_vector(blight_index_5.get_presence_query(my_query));
	//Presence_vector[0]==True if and only if the first kmer of the query sequence is in the graph
	
	//Compute the unique identifier of all the kmers of the query sequence
	vector<int64_t> hash_vector(blight_index_5.get_hashes_query(my_query));
	//hash_vector[0] indicate the identifier of the first kmer of the query 
	//A value of -1 mean that the kmer is not in the graph
	//Each identifier ID verify 0 <= ID < #kmers_in_DBG
	
	
	//We check the presence/absence of the query
	for(uint i(0);i<presence_vector.size();++i){
		if(presence_vector[i]){
			cout<<"Present ";
		}else{
			cout<<"Absent ";
		}
		
	}
	cout<<endl;
	
	
	//We check the presence/absence of the  bad query
	presence_vector=(blight_index_5.get_presence_query(my_query_bad));
	for(uint i(0);i<presence_vector.size();++i){
		if(presence_vector[i]){
			cout<<"Present ";
		}else{
			cout<<"Absent ";
		}
		
	}
	cout<<endl;
	
	
	//We print the identifier of the query
	for(uint i(0);i<hash_vector.size();++i){
		cout<<hash_vector[i]<<' ';
	}
	cout<<endl;
	
	
	//We print the identifier of the bad query
	hash_vector=blight_index_5.get_hashes_query(my_query_bad);
	for(uint i(0);i<hash_vector.size();++i){
		cout<<hash_vector[i]<<' ';
	}
	cout<<endl;
	
	

	//INDEX ON DISK
	
	
	
	//Write the index as a file on the disk
	blight_index_5.dump_disk("blight_index.gz");
	
	//Write the index as a file on the disk and deallocate the memory
	blight_index_5.dump_and_destroy("blight_index2.gz");
	
	//Load an index previously wrote on the disk
	kmer_Set_Light blight_index_6("blight_index.gz");
	
	
	//We print the identifier of the query performed on the loaded index
	hash_vector=blight_index_6.get_hashes_query(my_query);
	for(uint i(0);i<hash_vector.size();++i){
		cout<<hash_vector[i]<<' ';
	}
	cout<<endl;
	
	
	
	
	
	return 0;
}
