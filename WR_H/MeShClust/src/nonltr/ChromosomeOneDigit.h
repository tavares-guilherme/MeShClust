/*
 * ChromosomeOneDigit.h
 *
 *  Created on: Jul 31, 2012
 *      Author: Hani Zakaria Girgis, PhD - NCBI/NLM/NIH
 */

#ifndef CHROMOSOMEONEDIGIT_H_
#define CHROMOSOMEONEDIGIT_H_

#include <map>
#include "Chromosome.h"

namespace nonltr {
class ChromosomeOneDigit: public Chromosome {

private:
	/* Fields */
	map<char, int> * codes;	
	/*List of ints that represents the sequence of states*/
  	vector<int> *formatted_states;

	/* Methods */
	void help();
	void buildCodes();
	void encodeNucleotides();
	void printFormattedStates();

	void makeReverse();
	void makeComplement();
	void reverseSegments();

public:
	/* Methods */
	ChromosomeOneDigit();
	ChromosomeOneDigit(string);
	ChromosomeOneDigit(string, string);
	virtual ~ChromosomeOneDigit();
	virtual void finalize();
	vector<int> getFormattedStates();

	void makeR();
	void makeRC();
};
}

#endif /* CHROMOSOMEONEDIGIT_H_ */
