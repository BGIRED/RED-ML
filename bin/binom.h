/*
 * binom.h
 *
 *  Created on: 2014-10-26
 *      Author: Administrator
 */

#ifndef BINOM_H_
#define BINOM_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include "gzstream.h"

using namespace std;

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned long ulong;

void usage();
vector<string> split(string &);
double lngamm(int);
double lnfact(int);
double lnbico(int, int);
double binom(int, int, double);
void build_look(vector<double> &, uint, double, uint &);
void indexFinder(uint, uint, double, double, uint &);

void indexFinder(uint refNum, uint altNum, double falseRate, double biop, uint &index)
{
	uint depth(0),median(0);
	double mp(0.0);
	depth=refNum+altNum;
	if(refNum==0){
		index=0;
	}else if(refNum!=0){
		vector<double> table;
		build_look(table, depth, biop, median);
		mp = falseRate*altNum/refNum;
		for (index = median; index <= depth; ++index){
			if (table[index] < mp){
				break;
			}
		}
	}
}

void build_look(vector<double> &table, uint depth, double pvalue, uint &median){
	uint idx(0);
	double max(0.0);
	for (uint j(0); j <= depth; ++j){
		double p=binom(depth, j, pvalue);
		if(p>max){max=p;idx=j;}
		table.push_back(p);
	}
	median=idx;
}

/* IMPLEMENTATION */
/* ln gamma */
double lngamm(int z)
// Reference: "Lanczos, C. 'A precision approximation
// of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
// Translation of  Alan Miller's FORTRAN-implementation
// See http://lib.stat.cmu.edu/apstat/245
{
  double x = 0;
  x += 0.1659470187408462e-06/(z+7);
  x += 0.9934937113930748e-05/(z+6);
  x -= 0.1385710331296526    /(z+5);
  x += 12.50734324009056     /(z+4);
  x -= 176.6150291498386     /(z+3);
  x += 771.3234287757674     /(z+2);
  x -= 1259.139216722289     /(z+1);
  x += 676.5203681218835     /(z);
  x += 0.9999999999995183;
  return(log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5));
}

/* log factorial */
double lnfact(int n)
{
  if(n<=1) return(0);
  return(lngamm(n+1));
}

/* log binormal coefficient*/
double lnbico(int n, int k)
{
  return(lnfact(n)-lnfact(k)-lnfact(n-k));
}

double binom(int n, int k, double p){
	return exp(lnbico(n, k)) * pow(p, k) * pow(1-p, n-k);
}


#endif /* BINOM_H_ */
