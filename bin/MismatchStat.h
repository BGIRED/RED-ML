/*
 * MismatchStat.h
 *
 *  Created on: 2014-10-22
 *      Author: Administrator
 */

#ifndef MISMATCHSTAT_H_
#define MISMATCHSTAT_H_

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#include "SamCtrl.h"

using namespace std;

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned long long int ullint;

class parameter;
void usage();
void bamRead();
vector<string> split(string &);
inline string chr_check(const string &);
void callMismatch(string&,int&); // catch the mismatch information from the MD tag in bam file
void dump(ullint , ullint , ullint , ullint , ofstream &);

class parameter{
public:
	int minMapQ;
	string infile;
	string mode;
	string	outfile;
	int uniq;
    ullint maxRead;

	parameter(int _minMapQ, const string & _infile, const string & _mode, const string & _outfile, int _uniq, int _maxRead):
		minMapQ(_minMapQ), infile(_infile), mode(_mode), outfile(_outfile), uniq(_uniq), maxRead(_maxRead) {};
};


#endif /* MISMATCHSTAT_H_ */
