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
	
	
	
	//BLIGHT AS DICTIONARY/HASHMAP
	
	
	
	//We know the amont of distinct kmer present in the index
	cout<<"The index contains "<<blight_index_6.get_kmer_number()<<" distinct kmers"<<endl;
	
	//We allocate an  integer for each  kmer of the index and set them to 0
	uint abundance[blight_index_6.get_kmer_number()]={0};
	
	//We want to count the occurences of the indexed kmers in those dummy sequences
	string seq1("CTGATCGATCGTACGTAGCTGCTGATCGATCGTACGTACGTACGTCAGT");
	string seq2("CTGATCGATCGTACGTAGCTGCTGATCGATCGTA");
	string seq3("CTGATCGATCGTACGTAGCTGCTGATCGATC");
	string seq4("CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC");
	
	//We compute the indices of each kmer of the sequence
	hash_vector=blight_index_6.get_hashes_query(seq1);
	
	//Foreach computed indice
	for(uint i(0);i<hash_vector.size();++i){
		int indice(hash_vector[i]);
		
		//If the kmer is in the index its indice will be a position in the abundance vector
		if(indice!=-1){
			
			//We increment the associated counter
			abundance[indice]++;
		}
	}
	
	//We do so for the other sequences
	hash_vector=blight_index_6.get_hashes_query(seq2);
	for(uint i(0);i<hash_vector.size();++i){
		int indice(hash_vector[i]);
		if(indice!=-1){
			abundance[indice]++;
		}
	}
	
	hash_vector=blight_index_6.get_hashes_query(seq3);
	for(uint i(0);i<hash_vector.size();++i){
		int indice(hash_vector[i]);
		if(indice!=-1){
			abundance[indice]++;
		}
	}
	
	hash_vector=blight_index_6.get_hashes_query(seq4);
	for(uint i(0);i<hash_vector.size();++i){
		int indice(hash_vector[i]);
		if(indice!=-1){
			abundance[indice]++;
		}
	}
	
	//We take a look at the  computed abundances
	hash_vector=blight_index_6.get_hashes_query(my_query);
	for(uint i(0);i<hash_vector.size();++i){
		int indice(hash_vector[i]);
		if(indice!=-1){
			cout<<abundance[hash_vector[i]]<<' ';
		}
	}
	cout<<endl;
	
	
	
	//KMER INTERATOR
	
	
	
	//Another way to explore the index kmer is to use an iterator
	kmer_Set_Light_iterator it(&blight_index_6);
	
	vector<uint64_t> all_indices;
	do{
		//We can obtain a binary representation of the kmer as a integer
		kmer kmer_binary(it.get_kmer());
		cout<<"Kmer as an integer: "<<kmer_binary<<endl;
		
		//Or as a sequence
		string kmer_sequence(it.get_kmer_str());
		cout<<"Kmer in ASCII " <<kmer_sequence<<endl;
		
		//We compute the indice of each kmer and add them in a vector
		hash_vector=blight_index_6.get_hashes_query(kmer_sequence);
		all_indices.push_back(hash_vector[0]);
		
		//If iterator.next() return false it mean that there is no more kmer in the index
	}while(it.next());
	
	
	//We want to check that all the integers between 0 and number_kmer-1 are indice of an indexed kmer (Bijection)
	
	//We sort all the  stored indices
	sort(all_indices.begin(),all_indices.end());
	
	//It should be the list of integer [0,kmer_number[
	for(uint i(0);i<blight_index_6.get_kmer_number();++i){
		if(all_indices[i]!=i){
			cout<<"Noooo  the indices are not bijective with [O,number_kmer()[ !!!"<<endl;
			return 1;
		}else{
			cout<<i<<' ';
		}
	}
	cout<<"Everything went fine!"<<endl;
	
	return 0;
}
