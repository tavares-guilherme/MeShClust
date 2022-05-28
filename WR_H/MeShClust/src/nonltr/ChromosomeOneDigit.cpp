/*
 * ChromosomeOneDigit.cpp
 *
 *  Created on: Jul 31, 2012
 *      Author: Hani Zakaria Girgis, PhD at the NCB1/NLM/NIH
 * A	A
 * T	T
 * G	G
 * C	C
 * R	G or A
 * Y	T or C
 * M	A or C
 * K	G or T
 * S	G or C
 * W	A or T
 * H	A or C or T
 * B	G or T or C
 * V	G or C or A
 * D	G or T or A
 * N	G or T or A or C
 */
#include <iostream>
#include <map>
#include <list>

#include "Chromosome.h"
#include "ChromosomeOneDigit.h"
#include "../exception/InvalidInputException.h"

using namespace exception;

namespace nonltr {

ChromosomeOneDigit::ChromosomeOneDigit() :
		Chromosome() {
}

ChromosomeOneDigit::ChromosomeOneDigit(string fileName) :
		Chromosome(fileName) {
	help();
}

ChromosomeOneDigit::ChromosomeOneDigit(string seq, string info) :
		Chromosome(seq, info) {
	help();
}

void ChromosomeOneDigit::help() {
	// Build codes
	buildCodes();
	// Modify the sequence in the super class
	encodeNucleotides();
}

void ChromosomeOneDigit::finalize() {
	Chromosome::finalize();
	help();
}

void ChromosomeOneDigit::buildCodes() {
	// Make map
	codes = new map<char, int>();

	// Certain nucleotides
	codes->insert(map<char, char>::value_type('0', (int) 0));
	codes->insert(map<char, char>::value_type('1', (int) 1));
	codes->insert(map<char, char>::value_type('2', (int) 2));
	codes->insert(map<char, char>::value_type('3', (int) 3));
	codes->insert(map<char, char>::value_type('4', (int) 4));
	codes->insert(map<char, char>::value_type('5', (int) 5));
	codes->insert(map<char, char>::value_type('6', (int) 6));
	codes->insert(map<char, char>::value_type('7', (int) 7));
	codes->insert(map<char, char>::value_type('8', (int) 8));
	codes->insert(map<char, char>::value_type('9', (int) 9));

	// Common uncertain nucleotide
	// codes->insert(map<char, char>::value_type('N', (char) 4));

	/* Uncertain nucleotides
	codes->insert(map<char, char>::value_type('R', codes->at('G')));
	codes->insert(map<char, char>::value_type('Y', codes->at('C')));
	codes->insert(map<char, char>::value_type('M', codes->at('A')));
	codes->insert(map<char, char>::value_type('K', codes->at('T')));
	codes->insert(map<char, char>::value_type('S', codes->at('G')));
	codes->insert(map<char, char>::value_type('W', codes->at('T')));
	codes->insert(map<char, char>::value_type('H', codes->at('C')));
	codes->insert(map<char, char>::value_type('B', codes->at('T')));
	codes->insert(map<char, char>::value_type('V', codes->at('A')));
	codes->insert(map<char, char>::value_type('D', codes->at('T')));
	codes->insert(map<char, char>::value_type('N', codes->at('C')));
	codes->insert(map<char, char>::value_type('X', codes->at('G')));*/
}

ChromosomeOneDigit::~ChromosomeOneDigit() {
	codes->clear();
	delete codes;
}

/**
 * This method converts nucleotides in the segments to single digit codes
 */
/*	void ChromosomeOneDigit::encodeNucleotides() {
	cout << "Printing segment List" << endl;
	this->printSegmentList();
	cout << "segment printed" << endl;

	int aux;
	char state[3];
	int count;

	for (int s = 0; s < segment->size(); s++) {
		// Reset auxiliar variables
		count = 0;
		state[0] = '0';state[1] = '0'; state[2] = '0';
		
		int segStart = segment->at(s)->at(0);
		int segEnd = segment->at(s)->at(1);

		for (int i = segStart; i <= segEnd; i++) {			
			
			if (codes->count(base[i]) > 0) {
				base[i] = codes->at(base[i]);
			} else {
				cout << "reached here\n";
				string msg = "Invalid nucleotide: ";
				msg.append(1, base[i]);
				throw InvalidInputException(msg);
			}

			/* MyCode
			if (codes->count(base[i]) > 0) {
				state[count] = base[i];
				count++;
			} else if (base[i] == ','){
				aux = atoi(state);
				formatted_states.push_back(aux);
			}else{
				string msg = "Invalid nucleotide: ";
				msg.append(1, base[i]);
				throw InvalidInputException(msg);
			}
		}
  	}

	// Set last base of the sequence:
	//aux = atoi(state);
	//formatted_states.push_back(aux);

  // Digitize skipped segments
  int segNum = segment->size();
  if(segNum > 0){
    // The first interval - before the first segment
    int segStart = 0; 
    int segEnd = formatted_states.size(); 

    for (int s = 0; s <= segNum; s++) {      
		
		for (int i = segStart; i <= segEnd; i++) {
			char c = base[i];
			if(c != 'N'){
				if (codes->count(c) > 0) {
					base[i] = codes->at(c);
				} else {
					string msg = "Invalid nucleotide: ";
					msg.append(1, c);
					throw InvalidInputException(msg);
				}
			}
		}

		// The regular intervals between two segments
		if(s < segNum-1){
		segStart = segment->at(s)->at(1)+1;
		segEnd = segment->at(s+1)->at(0)-1;
		}
		// The last interval - after the last segment
		else if(s == segNum - 1){
		segStart = segment->at(s)->at(1)+1;
		segEnd = base.size()-1;
		} 
    } 
  }
}*/
void ChromosomeOneDigit::encodeNucleotides() {
	int seqLen = base.size();
	int aux;
	char state[4];
	int count = 0;
	
	formatted_states = new vector<int>;
	this->printFormattedStates();
	for (int i = 0; i < seqLen; i++) {
		
		if (codes->count(base[i]) > 0) {
			state[count] = base[i];
			count++;
		} else if (base[i] == ','){
			state[count] = '\0';
			aux = atoi(state);
			formatted_states->push_back(aux);
			count=0;
		}else{
			string msg = "Invalid nucleotide: ";
			msg.append(1, base[i]);
			throw InvalidInputException(msg);
		}
	}

	state[count] = '\0';
	aux = atoi(state);
	formatted_states->push_back(aux);


	this->printFormattedStates();

}

void ChromosomeOneDigit::printFormattedStates() {
	int l = formatted_states->size();

	cout << "State Sequences: = " << l << endl;
	for(int i = 0; i < l; i++){
		cout << formatted_states->at(i) << "\t";
	}
	cout << endl;
}

/**
 * Cannot be called on already finalized object.
 */
void ChromosomeOneDigit::makeR() {
	//cout << "Making reverse ..." << endl;
	makeReverse();
	reverseSegments();
}

/**
 * Cannot be called on already finalized object.
 */
void ChromosomeOneDigit::makeRC() {
	//cout << "Making reverse complement ..." << endl;
	makeComplement();
	makeReverse();
	reverseSegments();
}

void ChromosomeOneDigit::makeComplement() {
	map<char, int> complement;

	// Certain nucleotides
	complement.insert(map<char, int>::value_type((char) 0, (int) 4));
	complement.insert(map<char, int>::value_type((char) 1, (int) 3));
	complement.insert(map<char, int>::value_type((char) 2, (int) 2));
	complement.insert(map<char, int>::value_type((char) 3, (int) 1));
	complement.insert(map<char, int>::value_type((char) 4, (int) 0));

	// Unknown nucleotide
	complement.insert(map<char, char>::value_type('N', 'N'));
	// complement.insert(map<char, int>::value_type((char) 4, (int) 4));

	// Convert a sequence to its complement
	int seqLen = base.size();
	for (int i = 0; i < seqLen; i++) {
		if (complement.count(base[i]) > 0) {
			base[i] = complement.at(base[i]);
		} else {
			cerr << "Error: The digit " << (char) base[i];
			cerr << " does not represent a base." << endl;
			exit(2);
		}
	}
}

void ChromosomeOneDigit::makeReverse() {
	int last = base.size() - 1;

	// Last index to be switched
	int middle = base.size() / 2;

	for (int i = 0; i < middle; i++) {
		char temp = base[last - i];
		base[last - i] = base[i];
		base[i] = temp;
	}
}

void ChromosomeOneDigit::reverseSegments() {
	int segNum = segment->size();
	int lastBase = size() - 1;

	// Calculate the coordinate on the main strand
	for (int i = 0; i < segNum; i++) {
		vector<int> * seg = segment->at(i);

		int s = lastBase - seg->at(1);
		int e = lastBase - seg->at(0);
		seg->clear();
		seg->push_back(s);
		seg->push_back(e);
	}

	// Reverse the regions within the list
	int lastRegion = segNum - 1;
	int middle = segNum / 2;
	for (int i = 0; i < middle; i++) {
		vector<int> * temp = segment->at(lastRegion - i);
		(*segment)[lastRegion - i] = segment->at(i);
		(*segment)[i] = temp;
	}
}

}
