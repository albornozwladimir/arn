#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "main.h"
#include "utils.h"
#include "loader.h"
#include "options.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "constraints.h"
#include "traceback.h"
#include "subopt_traceback.h"
using namespace std;



void init_fold(string seq) {
	int len = seq.length();

	init_global_params(len);

	if (!encodeSequence(seq)) {
		free_fold(seq.length());
		exit(0);
	}
	
	create_tables(len);
	
	if (CONS_ENABLED) {
		init_constraints(constraintsFile.c_str(), len);
	}
}

void free_fold(int len) {
	if (CONS_ENABLED) 
		free_constraints(len);

	free_tables(len);
	free_global_params();
}

/**
 * Read the sequence out of the given filename and store it in seq
 *
 * @param filename A c string with the file to open
 * @param seq A C++ string object (passed by reference) to save to
 * @return SUCCESS or FAILURE
 */
int read_sequence_file(const char* filename, std::string& seq) {
	seq = "";

	ifstream fs;
	fs.open(filename, ios::in);
	if (fs == NULL)
		return FAILURE;

	string line;
	getline(fs, line);
	while(line.length() > 0) {
		// exclude lines starting with FASTA comment characters
		if(line[0] != ';' && line[0] != '>')
			seq += line;
		getline(fs, line);
	}

	fs.close();

    size_t loc;
    while((loc = seq.find(" ")) != string::npos)
        seq.erase(loc, 1);

	return SUCCESS;
}

int main(int argc, char** argv) {
	std::string seq;
	int energy;
	double t1;
	parse_options(argc, argv);

	if (read_sequence_file(seqfile.c_str(), seq) == FAILURE) {
		printf("Failed to open sequence file: %s.\n\n", seqfile.c_str());
		exit(-1);
	}
	
	// Read in thermodynamic parameters. Always use Turner99 data (for now)
	readThermodynamicParameters("Turner99",false);
	printRunConfiguration(seq);
	init_fold(seq);
	printf("\nComputing minimum free energy structure...\n");
	fflush(stdout);
	energy = calculate(seq.length(), nThreads);
	printf("Done.\n\n");
	printf("Results:\n");
	printf("- Minimum Free Energy: %12.4f kcal/mol\n", energy/100.00);
	printf("- MFE runtime: %9.6f seconds\n", t1);
	
	printf("\n");
	print_sequence(seq.length());
	print_structure(seq.length());
	if (CONS_ENABLED)
		print_constraints(seq.length());

	// release the malloc'd arrays
	free_fold(seq.length());

	return EXIT_SUCCESS;
}
