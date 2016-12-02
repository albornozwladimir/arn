#ifndef _MAIN_H
#define _MAIN_H

#include <string>
#include <vector>
#include "constants.h"

using namespace std;

//void traverse_wholeData(int[], int, int);

enum GTFOLD_FLAGS { SUCCESS = 0, FAILURE, ERR_OPEN_FILE, NO_CONS_FOUND};

GTFOLD_FLAGS initialize_constraints(int*** fbp, int*** pbp, int& numpConstraints, int& numfConstraints, const char* constr_file);

GTFOLD_FLAGS handle_IUPAC_code(const std::string& s, const int bases);
void limit_contact_distance(int lCD, int length);
void force_noncanonical_basepair(const char* nc_basepairs, int length);

bool is_valid_base(char c)
{	
	return ( (c-'A' == 0) || (c-'a' == 0) || 
			 (c-'C' == 0) || (c-'c' == 0) ||
			 (c-'G' == 0) || (c-'g' == 0) ||
			 (c-'U' == 0) || (c-'u' == 0));
}

void trim_spaces(std::string& str)
{
	// Trim Both leading and trailing spaces
	size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
	size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( string::npos == startpos ) || ( string::npos == endpos))
	{
		str = "";
	}
	else
		str = str.substr( startpos, endpos-startpos+1 );

}

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ")
{
	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (std::string::npos != pos || std::string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}
	
	


#endif
