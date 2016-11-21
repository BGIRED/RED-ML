/*
 * MismatchStat.cpp
 *
 *  Created on: 2014-10-22
 *      Author: Administrator
 */
#include"MismatchStat.h"

#define PROGRAM_COMPILE_DATE __DATE__
#define PROGRAM_COMPILE_TIME __TIME__
parameter para(0,"","bam","",0,1000000000); // minMapQ,minBaseQ,baseQshift,infile,mode,ref,outfile,uniq,discard

int main(int argc, char** argv) {
	for (int i = 1; i < argc; i++) {
	    if(argv[i][0]=='-'){
	        switch (argv[i][1]) {
	        	case 'i': para.infile = argv[++i];break;
	        	case 'm': para.mode= argv[++i];break;
	        	case 'q': para.minMapQ = atoi(argv[++i]); break;
	        	case 'o': para.outfile = argv[++i];break;
                case 'x': para.maxRead = atoi(argv[++i]);break;
	        	case 'u': para.uniq = 1;break;
	            case 'h': usage();exit(0);
	            default : usage();exit(0);
	        }
	    }
	}
	if(para.infile.empty()) {usage();exit(0);}
	if(para.outfile.empty()) {usage();exit(0);}
	time_t timep;
	time(&timep);
	cout<<"Begin at: "<<ctime(&timep);
	bamRead();
	time(&timep);
	cout<<"End at: "<<ctime(&timep);
	return 0;
}

void bamRead()
{
		ofstream STAT;
		STAT.open(para.outfile.c_str());
		if ( ! STAT.good()) {
		    cerr << "ERROR: Opening file `" << para.outfile << "' failed.\n";
		    exit(0);
		}
		SamCtrl SAM;
		string fileMode;
		if(para.mode.compare("sam")==0 || para.mode.compare("s")==0) fileMode="r";
			else fileMode="rb";

		if(SAM.open(para.infile.c_str(),fileMode.c_str()) ) {
			cout<<"Reading " << para.infile << " ..." << endl;
		}else{
			cerr << "ERROR: Opening file `" << para.infile << "' failed.\n";
			exit(0);
		}

		string line;
		vector<string> bin;
		int mapQ(0);
		string cigar,chr,last_chr("");
		ullint lineNum(0);
		ullint uniqNum(0);
		ullint totMism(0);
		ullint totBase(0);

		while(SAM.readline(line)!= -1)
		{
			bin=split(line);
			mapQ=atoi(bin[4].c_str());
			if(mapQ <para.minMapQ ) continue;
			if(bin.size()<11 || bin[5] == "*" || bin[5].find('M')==string::npos) continue;
			int bestHit(0);
			string mds,seq;
            int flag=atoi(bin[1].c_str()); // 2015/1/1 by Donby
            if(flag & 0x0400) continue; // 2015/1/1 by Donby, remove PCR dup
			for (int u=11; u<bin.size(); u++){
				if(bin[u].find("NH:", 0) == 0 || bin[u].find("H0:", 0) == 0 || bin[u].find("X0:", 0) == 0 || bin[u].find("IH:", 0) == 0){ // using the unique mapped reads, 2015/1/7 by Donby
					bestHit = atoi(bin[u].substr(5, bin[u].size()-5).c_str());
				}else if(bin[u].find("MD:",0) == 0){
					mds = bin[u].substr(5, bin[u].size()-5);
				}
			}
			lineNum++;
			if(bestHit!=1 && para.uniq) continue; // 2015/1/7 by Donby
			if(para.uniq && bestHit==1) uniqNum++;
            if(para.uniq)
            {
                if(uniqNum>para.maxRead) // ensure to use enough uniq mapping reads when para.uniq is set. 2015/04/08 by Donby
                {
                    lineNum--;uniqNum--;
                    dump(uniqNum,lineNum,totBase,totMism,STAT);
                    exit(0);
                }
            }else{
                if(lineNum>para.maxRead) // use all the reads. 2015/04/08 by Donby
                {
                    lineNum--;
                    dump(uniqNum,lineNum,totBase,totMism,STAT);
                    exit(0);
                }
            }

			chr=chr_check(bin[2]);
			if(chr!=last_chr){
				cout << "Program traverse to: " << chr << endl;
				last_chr=chr;
			}

			int mism(0); // number of mismatches on a read
			seq=bin[9];
			int qlen=seq.size();
			callMismatch(mds,mism);
			totMism+=mism;
			totBase+=qlen;
//			cout<<chr<<"\t"<<bin[3]<<"\t"<<mism<<"\t"<<qlen<<"\t"<<totBase<<endl;
		}
        dump(uniqNum,lineNum,totBase,totMism,STAT);
		time_t timep;
		time(&timep);
		cout<<"End at: "<<ctime(&timep);
		STAT.close();
		SAM.close();
}

void dump(ullint uniqNum, ullint lineNum, ullint totBase, ullint totMism, ofstream &STAT)
{
    double mismRate;
    mismRate=static_cast<double>(totMism)/static_cast<double>(totBase);
    STAT<<"Mismatch bases    : "<<totMism<<endl;
    STAT<<"Total bases       : "<<totBase<<endl;
    STAT<<"Mismatch rate     : "<<mismRate<<endl;
    STAT<<"Unique read number: "<<uniqNum<<endl;
    STAT<<"Mapped read number: "<<lineNum<<endl;
}

void callMismatch(string &mds, int &mism)
{
	for (int m=0; m<mds.size(); m++)
	{
		if (mds[m] == '^')
		{
			for ( m++; m<mds.size(); m++)
			{
				if (isdigit(mds[m])) break;
			}
		}else{
			if (isalpha(mds[m])) mism++;
		}
	}
}

vector<string> split(string& line)
{
	vector<string> word;
	stringstream IN(line);
	string tmp;
	while (!IN.eof())
	{
		IN >> tmp;
		//cout<<tmp<<endl;
		word.push_back(tmp);
	}
	return word;
}

inline string chr_check(const string &bchr) {
	string chr;
	if (bchr.size()<4 || (bchr.size()>3 && bchr.substr(0,3)!="chr")){
		chr="chr";
		chr+=bchr;
	}else{
		chr=bchr;
	}
	return chr;
}

void usage(){
	cerr << endl;
	cerr << "BGI2294  "<<PROGRAM_COMPILE_DATE<<' '<<PROGRAM_COMPILE_TIME<<endl;
	cerr << endl;
	cerr << "Usage:  " <<"MismatchStat -i in.bam -o stat.txt -u" << endl;
	cerr << '\t' << "-i [str] input sorted alignment file, required" << endl;
	cerr << '\t' << "-m [str] format of input file, 'sam' for sam, 'bam' for bam, default[bam]" << endl;
	cerr << '\t' << "-q [int] minimal mapping quality, default[0]" << endl;
    cerr << '\t' << "-x [int] using the first int reads to calculate, default[1000000000]" << endl;
	cerr << '\t' << "-o [str] output file, required" << endl;
	cerr << '\t' << "-u use unique mapping reads to calculate mismatch rate, default[off]" << endl;
	cerr << '\t' << "-h print this usage" << endl;
	cerr << endl;
	return;
}




