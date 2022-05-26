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

	/* Methods */
	void help();
	void buildCodes();
	void encodeNucleotides();

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

	void makeR();
	void makeRC();
};
}

#endif /* CHROMOSOMEONEDIGIT_H_ */
