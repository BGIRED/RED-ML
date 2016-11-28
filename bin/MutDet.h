/*
 * pileup.h
 *
 *  Created on: 2014-9-20
 *      Author: Donby
 *     Contact: liudongbing@genomics.cn
 */

#ifndef PILEUP_H_
#define PILEUP_H_

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <sstream>
#include<map>
#include<cmath>
#include<algorithm>
#include <cstdlib>
#include "SamCtrl.h"
#include "gzstream.h"
#include"binom.h"
#include"fisher.h"

using namespace std;

typedef unsigned int uint;
typedef unsigned char uchar;

class parameter;
class cell;
class pileup;
void usage();
void bamRead();
void parseFa();
void dump(string,ogzstream &, pileup &);
void initial(string); // initial the pile class
vector<string> split(string &);
inline string chr_check(const string &);
void callMismatch(string&,int&); // catch the mismatch information from the MD tag in bam file
void readPileup(string &,int&,string&,string&,string&, string &, string &, string &, string &, pileup &, pileup &); // pileup a read
void display(string &,int &,int &, const int &, const string &, const string &, string &, string &, string &, string &, pileup &, pileup &); // store the detailed information for each site in the pile class
inline string inTstring(int &);
void getP();
void siteCalc(string, int, char, cell &, ogzstream & );
void maxMapValue(map<char,vector<uint> > &, char&);
void StringSplit(std::string s, char splitchar, std::vector<std::string>& vec);
void strtoivec(vector<string> &, vector<int> &);

map<string,string> refseq;
map<string,int> seqlen;
map<char,vector<uint> > alleleInfo;

class parameter{
public:
	int minMapQ;
	int minBaseQ;
	int minDepth;
	int minAlt;
	int rstart;
	int rend;
	int baseQshift;
	string infile;
	string mode;
	string ref;
	string	outfile;
	int uniq;
	double pvalue;
	double evalue;
	string statFile;
    int emitAll; // 2015/1/8 added by Donby

	parameter(int _minMapQ, int _minBaseQ, int _minDepth, int _minAlt, int _rstart, int _rend, int _baseQshift, const string & _infile, const string & _mode, const string & _ref, const string & _outfile, int _uniq, double _pvalue, double _evalue, string _statFile, int _emitAll):
		minMapQ(_minMapQ), minBaseQ(_minBaseQ), minDepth(_minDepth), minAlt(_minAlt), rstart(_rstart), rend(_rend), baseQshift(_baseQshift), infile(_infile), mode(_mode), ref(_ref), outfile(_outfile), uniq(_uniq), pvalue(_pvalue), evalue(_evalue), statFile(_statFile), emitAll(_emitAll) {};
};

// store the depth, bases, map quals, base quals, read position, mismatch number, read length information
class cell{
public:
	int depth;
	string bases;
	string mapQs;
	string baseQs;
	string qposs;
	string misms;
	string qlens;

	cell(): depth(0), bases(""), mapQs(""), baseQs(""), qposs(""), misms(""), qlens(""){}

	void operator=(const cell &source){
		depth = source.depth;
		bases = source.bases;
		mapQs = source.mapQs;
		baseQs = source.baseQs;
		qposs = source.qposs;
		misms = source.misms;
		qlens = source.qlens;
	}

	void clear()
	{
		depth=0;
		bases.clear();
		mapQs.clear();
		baseQs.clear();
		qposs.clear();
		misms.clear();
		qlens.clear();
	}
};

class pileup{
public:
	int start;
	int end;
	int size;
	int rsize;
	vector<cell> record;

	pileup(int s, int e):start(s), end(e), size(e - s + 1), rsize(0){record.resize(size);}

	void operator=(const pileup &source){
		if (&source != this){
			start=source.start;
			end=source.end;
			size=source.size;
			for (size_t i(0); i < source.record.size(); ++i)
				record[i] = source.record[i];
		}
	}
	inline void clear();
	inline void clear(const int);
	inline void resize(const int, const int);
	inline int copy(const pileup &, const int, const int);
	inline int icopy(const int, const int);
};

inline void pileup::clear(){
	for (size_t i(0); i < record.size(); ++i)
		record[i].clear();
	rsize = 0;
}

inline void pileup::clear(const int s){
	if (s >= rsize)
		return;

	for (int i(s); i < rsize; ++i)
		record[i].clear();

	rsize = s;
}


inline void pileup::resize(const int s, const int e){
	start = s;
	end = e;
	size = e - s + 1;
	rsize = 0;
	record.resize(size);
}

inline int pileup::copy(const pileup &source, const int s, const int e){
	int i(0), j(s);
	for(; j < e && i < size; ++j, ++i)
		record[i] = source.record[j];

	if (rsize < i)
		rsize = i;

	return j;
}

inline int pileup::icopy(const int s, const int e){
	int i(0);
	for (int j(s); j < e && j < rsize; ++j, ++i)
		record[i] = record[j];

	(*this).clear(e - s);
	rsize = i;

	return i;
}
#endif /* PILEUP_H_ */
