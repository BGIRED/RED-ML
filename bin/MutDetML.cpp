/*
 * pileup.cpp
 *
 *  Created on: 2014-9-20
 *      Author: Donby
       Contact: liudongbing@genomics.cn
 */
#include"MutDet.h"

#define PROGRAM_COMPILE_DATE __DATE__
#define PROGRAM_COMPILE_TIME __TIME__
#define BUFFER_SIZE  1000000
parameter para(20,20,4,2,5,10,33,"","bam","","",0,0.01,0.01,"",0); // minMapQ,minBaseQ,baseQshift,infile,mode,ref,outfile,uniq,discard

int main(int argc, char** argv) {
	for (int i = 1; i < argc; i++) {
	    if(argv[i][0]=='-'){
	        switch (argv[i][1]) {
	        	case 'i': para.infile = argv[++i];break;
	        	case 'm': para.mode= argv[++i];break;
	        	case 'r': para.ref = argv[++i];break;
	        	case 'q': para.minMapQ = atoi(argv[++i]); break;
	        	case 'Q': para.minBaseQ = atoi(argv[++i]);break;
	        	case 's': para.baseQshift = atoi(argv[++i]);break;
	        	case 'd': para.minDepth = atoi(argv[++i]);break;
	        	case 'a': para.minAlt = atoi(argv[++i]);break;
	        	case 'l': para.rstart = atoi(argv[++i]);break;
	        	case 'n': para.rend = atoi(argv[++i]);break;
	        	case 'e': para.evalue = strtod(argv[++i],NULL);break;
	        	case 'v': para.statFile = argv[++i];break;
	        	case 'o': para.outfile = argv[++i];break;
	        	case 'u': para.uniq = 1;break;
                case 't': para.emitAll = 1;break; // Added by Donby 2015/1/8
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
    cout << "Parameter setting:" << endl;
    cout << "\tInput file:\t" << para.infile << endl;
    cout << "\tOutput file:\t" << para.outfile << endl;
    cout << "\tMinimal base quality:\t" << para.minBaseQ << endl;
    cout << "\tMinimal mapping quality:\t" << para.minMapQ << endl;
    cout << "\tMinimal depth required:\t" << para.minDepth << endl;
    cout << "\tMinimal alt support:\t" << para.minAlt << endl;
    cout << "\tRead start offset:\t" << para.rstart << endl;
    cout << "\tRead end offset:\t" << para.rend << endl;
    cout << "\tFalse discovery rate:\t" << para.evalue << endl;
    if ( para.uniq ) cout<< "\tUnique mapping flag is set on" << endl;
    cout << "\tBackgroud mismatch rate file:\t" << para.statFile << endl;
    para.minBaseQ+=para.baseQshift;
	getP();
    if(para.emitAll){cout << "\tEmit all the sites" << endl;} // Added by Donby 2015/1/8
    if(!para.ref.empty()){parseFa();}
	bamRead();
	time(&timep);
	cout<<"End at: "<<ctime(&timep);
	return 0;
}

void bamRead()
{
		ogzstream PILEUP;
		PILEUP.open(para.outfile.c_str());
		if ( ! PILEUP.good()) {
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
		int pos;
		vector<string> bin;
		int mapQ(0);
		string cigar,seq,qual,chr;
		int lineNum(0);
		string last_chr("");

		int interval=10000000;
		int traverse(interval);
		int chrNum(0);
		pileup seg(1,BUFFER_SIZE);
		pileup buf(BUFFER_SIZE+1,2*BUFFER_SIZE);
        PILEUP<<"Chr\tPos\tDepth\tRef\tSup\tMid\tEnd\tPlus\tMinus\tMism\tNoMism\tMeanBQ\tAlt\tSup\tMid\tEnd\tPlus\tMinus\tMism\tNoMism\tMeanBQ\tIndex\trposPvalue\tsbPvalue\tmmPvalue"<<endl;

		while(SAM.readline(line)!= -1)
		{
			bin=split(line);
			mapQ=atoi(bin[4].c_str());
			if(mapQ <para.minMapQ ) continue;
			if(bin.size()<11 || bin[5] == "*" || bin[5].find('M')==string::npos) continue;
			int bestHit(0);
			string mds;
            int flag=atoi(bin[1].c_str()); // 2014/12/29 by Donby
            if(flag & 0x0400) continue; // 2014/12/29 by Donby, remove PCR dup
			for (int u=11; u<bin.size(); u++){
				if(bin[u].find("NH:", 0) == 0 || bin[u].find("H0:", 0) == 0 || bin[u].find("X0:", 0) == 0 || bin[u].find("IH:", 0) == 0){ // using the unique mapped reads, 2015/01/07 by Donby
					bestHit = atoi(bin[u].substr(5, bin[u].size()-5).c_str());
				}else if(bin[u].find("MD:",0) == 0){
					mds = bin[u].substr(5, bin[u].size()-5);
				}
			}
			if(bestHit!=1 && para.uniq) continue; // 2015/1/7 by Donby, remove ununique mapped reads
			if(seq.size()!=qual.size()) {
				cerr << "Read length do not equal Qual length at line: " << line << endl;
				continue;
			}

			chr=chr_check(bin[2]);
			pos=atoi(bin[3].c_str());
			if(chr!=last_chr)
			{
				chrNum++;
				if(!para.ref.empty()){
					if(refseq.find(chr)==refseq.end())
                    {
                        cerr << "Can not find " << chr << " in reference file" << endl;
                        continue;
                    }
				}
				if(chrNum!=1){
					dump(last_chr,PILEUP,seg);
					seg=buf;
					dump(last_chr,PILEUP,seg);
					seg.clear();
					buf.clear();
					traverse=interval;
				}
				int newStart=pos;
				int newEnd=pos+BUFFER_SIZE-1;
				seg.resize(newStart,newEnd);
				newStart=pos+BUFFER_SIZE;
				newEnd=pos+2*BUFFER_SIZE-1;
				buf.resize(newStart,newEnd);
				last_chr=chr;
				cout << "Processing " << chr << endl;
			}

			while(pos>=traverse)
			{
				cout << "Program traverse to: " << chr << ":" << traverse << endl;
				traverse+=interval;
			}
			while(pos-seg.start>=BUFFER_SIZE)
			{
				dump(last_chr,PILEUP,seg);
				seg=buf;
				buf.clear();
				buf.resize(seg.end+1,seg.end+BUFFER_SIZE);
			}

			cigar=bin[5];
			seq=bin[9];
			qual=bin[10];
//            int flag=atoi(bin[1].c_str());
			int mism(0); // number of mismatches on a read
			callMismatch(mds,mism);
			int qlen=seq.size();
			string smism=inTstring(mism);
			string smapQ=inTstring(mapQ);
			string sqlen=inTstring(qlen);
            string strand="+";
            if(flag & 0x0010) strand="-";
			readPileup(bin[0],pos,seq,qual,cigar,smapQ,smism,sqlen,strand,seg,buf);
		}

		dump(last_chr,PILEUP,seg);
		seg=buf;
		buf.clear();
		dump(last_chr,PILEUP,seg);
		seg.clear();
		buf.clear();
		time_t timep;
		time(&timep);
		cout<<"End at: "<<ctime(&timep);
		PILEUP.close();
		SAM.close();
}

void readPileup(string &readID, int &pos, string &seq, string &qual, string &cigar, string &smapQ, string &smism, string &sqlen, string &strand, pileup &seg, pileup &buf){
	int len = 0; // the length of digit of each cigar part
	int beg = 0; // the begin position on cigar
	int end; // the end position on cigar
	int qpos = 1; //the position on read
	int rpos = pos; // the index for object pile
	for (end=0; end<cigar.size(); end++){
		switch (cigar[end]){
			case 'M':
			{
				len= atoi(cigar.substr(beg, end - beg).c_str());
				display(readID,rpos,qpos,len,seq,qual,smapQ,smism,sqlen,strand,seg,buf);
				beg=end+1;
				rpos+=len;
				qpos+=len;
				break;
			}
			case 'I':
			{
				len= atoi(cigar.substr(beg, end - beg).c_str());
				beg=end+1;
				qpos+=len;
				break;
			}
			case 'D':
			{
				len= atoi(cigar.substr(beg, end - beg).c_str());
				beg=end+1;
				rpos+=len;
				break;
			}
			case 'N':
			{
				len= atoi(cigar.substr(beg, end - beg).c_str());
				beg=end+1;
				rpos+=len;
				break;
			}
			case 'S':
			{
				len= atoi(cigar.substr(beg, end - beg).c_str());
				beg=end+1;
				qpos+=len;
				break;
			}
			case 'H':
			{
				len= atoi(cigar.substr(beg, end - beg).c_str());
				beg=end+1;
				break;
			}
			case 'P':
			{
				len= atoi(cigar.substr(beg, end - beg).c_str());
				beg=end+1;
				break;
			}
			default : break;
		}
	}
}

void display(string &readID,int &rpos,int &qpos, const int &len, const string &seq, const string &qual, string &smapQ, string &smism, string &sqlen, string &strand, pileup &seg, pileup &buf){
	for (int p=0; p<len; p++){
		int rpidx=rpos+p;
		int qpidx=qpos+p;
		int ridx=rpidx-1;
		int qidx=qpidx-1;
		int sidx=rpidx-seg.start;
        int bidx=rpidx-buf.start;
		int Qual=int(qual[qidx]);
		if(Qual<para.minBaseQ) continue;
		string sbase="";
        if(strand=="-")
        {
            sbase=tolower(seq[qidx]);
        }else{
            sbase=toupper(seq[qidx]);
        }
        if(sbase=="N") continue;
		string squal=qual.substr(qidx,1);
        string sqpidx=inTstring(qpidx);
		if(rpidx>=seg.start && rpidx<=seg.end)
		{
			seg.record[sidx].depth++;
            seg.record[sidx].bases += sbase;
            seg.record[sidx].baseQs += squal;
            if(seg.record[sidx].depth==1)
            {
		    	seg.record[sidx].mapQs += smapQ;
			    seg.record[sidx].qposs += sqpidx;
    			seg.record[sidx].misms += smism;
	    		seg.record[sidx].qlens += sqlen;
            }else{
                seg.record[sidx].mapQs += ",";
                seg.record[sidx].mapQs += smapQ;
                seg.record[sidx].qposs += ",";
                seg.record[sidx].qposs += sqpidx;
                seg.record[sidx].misms += ",";
                seg.record[sidx].misms += smism;
                seg.record[sidx].qlens += ",";
                seg.record[sidx].qlens += sqlen;
            }
		}else if(rpidx>=buf.start && rpidx<=buf.end)
		{
            buf.record[bidx].depth++;
            buf.record[bidx].bases += sbase;
            buf.record[bidx].baseQs += squal;
            if(buf.record[bidx].depth==1)
            {
                 buf.record[bidx].mapQs += smapQ;
                 buf.record[bidx].qposs += sqpidx;
                 buf.record[bidx].misms += smism;
                 buf.record[bidx].qlens += sqlen;
            }else{
                buf.record[bidx].mapQs += ",";
                buf.record[bidx].mapQs += smapQ;
                buf.record[bidx].qposs += ",";
                buf.record[bidx].qposs += sqpidx;
                buf.record[bidx].misms += ",";
                buf.record[bidx].misms += smism;
                buf.record[bidx].qlens += ",";
                buf.record[bidx].qlens += sqlen;
            }
		}else{
			cerr<<"Read spaned longer than "<<BUFFER_SIZE<<" at "<<readID<<"\t"<<rpos<<endl;
		}
	}
}

void dump(string chr, ogzstream & OP, pileup &seg)
{
	for(size_t i=0;i!=seg.record.size();i++)
	{
		size_t j=i+seg.start;
		if (seg.record[i].depth == 0) continue;
		if (j>seqlen[chr]) continue;
		char ref=toupper(refseq[chr][j-1]);
//		OP<<chr<<"\t"<<j<<"\t"<<ref<<"\t"<<seg.record[i].depth<<"\t"<<seg.record[i].bases<<"\t"<<seg.record[i].baseQs<<"\t"<<seg.record[i].mapQs<<"\t"<<seg.record[i].qposs<<"\t"<<seg.record[i].misms<<"\t"<<seg.record[i].qlens<<endl;
		cell siteInfo;
		siteInfo=seg.record[i];
		siteCalc(chr,j,ref,siteInfo,OP);
	}
}

void siteCalc(string chr, int pos, char ref, cell &siteInfo, ogzstream & OP)
{
	alleleInfo.clear();
	vector<uint> siteStat(8,0);//sup,mid,end,plus,minus,mism,no-mism
	vector<string> qposs;
	vector<string> qlens;
	vector<string> misms;
	StringSplit(siteInfo.qposs, ',', qposs);
	StringSplit(siteInfo.qlens, ',', qlens);
	StringSplit(siteInfo.misms, ',', misms);
	vector<int> qposi;
	vector<int> qleni;
	vector<int> mismi;
	strtoivec(qposs,qposi);
	strtoivec(qlens,qleni);
	strtoivec(misms,mismi);
	for (int i=0,j=0; i<siteInfo.depth; i++,j+=2)
	{
		char sallele(siteInfo.bases[i]);
		char allele(toupper(sallele));
		int alleleq=int(siteInfo.baseQs[i])-para.baseQshift;
		if(alleleInfo.find(allele)==alleleInfo.end())
		{
			alleleInfo.insert(make_pair(allele,siteStat));
		}
		alleleInfo[allele][0]++; // allele depth
		alleleInfo[allele][7]+=alleleq; // sum of allele qual (qualified)
		int qpos,qlen,mism,end;
		qpos=qposi[i];
		qlen=qleni[i];
		mism=mismi[i];
//		cout<<chr<<"\t"<<pos<<"\t"<<ref<<"\t"<<j<<"\t"<<sqpos<<"\t"<<sqlen<<"\t"<<smism<<endl;
		if(isupper(sallele))
		{
			alleleInfo[allele][3]++;
		}else{
			alleleInfo[allele][4]++;
			qpos=qlen-qpos+1;
		}
		if(allele==ref)
		{
			if(mism==0)
			{
				alleleInfo[allele][6]++;
			}else{
				alleleInfo[allele][5]++;
			}
		}else{
			if(mism<=1)
			{
				alleleInfo[allele][6]++;
			}else{
				alleleInfo[allele][5]++;
			}
		}
		if(qpos<=qlen/2)
		{
			if(qpos<=para.rstart)
			{
				alleleInfo[allele][2]++;
			}else{
				alleleInfo[allele][1]++;
			}
		}else{
			if((qlen-qpos+1)<=para.rend)
			{
				alleleInfo[allele][2]++;
			}else{
				alleleInfo[allele][1]++;
			}
		}
	}

	vector<uint> refStat(8,0);
	for (map<char,vector<uint> >::iterator ita=alleleInfo.begin(); ita!=alleleInfo.end(); ita++)
	{
		if(ita->first==ref){
			refStat=ita->second;
			alleleInfo.erase(ita);
			break;
		}
	}
	char alt;
	vector<uint> altStat(8,0);
	if(alleleInfo.empty()){
        if(para.emitAll){
			double meanRefQ(0.0);
			meanRefQ=double(refStat[7])/double(refStat[0]);
            OP<<chr<<"\t"<<pos<<"\t"<<refStat[0]<<"\t"<<ref<<"\t"<<refStat[0]<<"\t"<<refStat[1]<<"\t"<<refStat[2]<<"\t"<<refStat[3]<<"\t"<<refStat[4]<<"\t"<<refStat[5]<<"\t"<<refStat[6]<<"\t"<<meanRefQ<<"\tN\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\t1"<<endl; // 2015/1/8 added by Donby
        }
	}else{
		maxMapValue(alleleInfo,alt);
		altStat=alleleInfo[alt];
		double readPosP(1.0),readPosL(1.0),readPosR(1.0),strandP(1.0),strandL(1.0),strandR(1.0),mismP(1.0),mismL(1.0),mismR(1.0),meanRefQ(0.0),meanAltQ(0.0); // 2015/1/8 modified from 0.0 to 1.0 by Donby, 2015/12/01 add meanRefQ and meanAltQ
		uint index(0),refNum(0),altNum(0),depth(0);
		refNum=refStat[0];
		altNum=altStat[0];
		depth=refNum+altNum;
		meanRefQ=double(refStat[7])/double(refNum);
		meanAltQ=double(altStat[7])/double(altNum);
		if(altNum>=para.minAlt && depth>=para.minDepth)
		{
            if(refNum==0)
            {
                index=0;
            }else{
    			indexFinder(refNum, altNum, para.evalue, para.pvalue, index);
            }
			kt_fisher_exact(refStat[1], altStat[1], refStat[2], altStat[2] , &readPosL, &readPosR, &readPosP);
			kt_fisher_exact(refStat[3], altStat[3], refStat[4], altStat[4] , &strandL, &strandR, &strandP);
			kt_fisher_exact(refStat[5], altStat[5], refStat[6], altStat[6] , &mismL, &mismR, &mismP);
			OP<<chr<<"\t"<<pos<<"\t"<<depth<<"\t"<<ref<<"\t"<<refNum<<"\t"<<refStat[1]<<"\t"<<refStat[2]<<"\t"<<refStat[3]<<"\t"<<refStat[4]<<"\t"<<refStat[5]<<"\t"<<refStat[6]<<"\t"<<meanRefQ<<"\t"<<alt<<"\t"<<altNum<<"\t"<<altStat[1]<<"\t"<<altStat[2]<<"\t"<<altStat[3]<<"\t"<<altStat[4]<<"\t"<<altStat[5]<<"\t"<<altStat[6]<<"\t"<<meanAltQ<<"\t"<<index<<"\t"<<readPosP<<"\t"<<strandP<<"\t"<<mismP<<endl; // 2015/12/01 add meanRefQ and meanAltQ
		}else{
            if(para.emitAll){
                OP<<chr<<"\t"<<pos<<"\t"<<depth<<"\t"<<ref<<"\t"<<refNum<<"\t"<<refStat[1]<<"\t"<<refStat[2]<<"\t"<<refStat[3]<<"\t"<<refStat[4]<<"\t"<<refStat[5]<<"\t"<<refStat[6]<<"\t"<<meanRefQ<<"\t"<<alt<<"\t"<<altNum<<"\t"<<altStat[1]<<"\t"<<altStat[2]<<"\t"<<altStat[3]<<"\t"<<altStat[4]<<"\t"<<altStat[5]<<"\t"<<altStat[6]<<"\t"<<meanAltQ<<"\t"<<index<<"\t"<<readPosP<<"\t"<<strandP<<"\t"<<mismP<<endl; // 2015/1/8 added by Donby; 2015/12/01 add meanRefQ and meanAltQ
            }
        }
	}
}

void maxMapValue(map<char,vector<uint> > &alleleInfo, char &maxAlt)
{
	uint max(0),flag(0);
	for (map<char,vector<uint> >::iterator ita=alleleInfo.begin(); ita!=alleleInfo.end(); ita++)
	{
		flag++;
		if(flag==1)
		{
			max=ita->second[0];
			maxAlt=ita->first;
		}else{
			if(ita->second[0]>max)
			{
				max=ita->second[0];
				maxAlt=ita->first;
			}
		}
	}
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

void getP()
{
    ifstream STAT(para.statFile.c_str());
    if(!STAT.good()){
        cerr << "ERROR: Opening file `" << para.statFile << endl;
        exit(0);
    }
    string line,pv;
    int flag(0);
    while (getline(STAT,line))
    {
        vector<string> vline;
        vline=split(line);
        if(vline[0]=="Mismatch" && vline[1]=="rate")
        {
            para.pvalue=strtod(vline[3].c_str(),NULL);
            cout << "\tBackground mismatch rate:\t" << para.pvalue << endl;
            flag++;
        }
    }
    if(flag==0)
    {
        cerr << "Couldn't find Mismatch rate in file: " << para.statFile << endl;
        exit(0);
    }else if(flag!=1){
        cerr << "Mismatch rate finding error: " << para.statFile <<endl;
        exit(0);
    }
}

void parseFa()
{
		igzstream FA(para.ref.c_str());
		if ( ! FA.good()) {
		    cerr << "ERROR: Opening file `" << para.ref << endl;
		    exit(0);
		}else{
			cout << "Reading " << para.ref << " ..." << endl;
		}
		string instr,chr,seq;
		vector<string> bin;
		int n=0;
		while (getline(FA,instr)) {
			if (instr[0]=='>'){
				n++;
				if (n>1)
				{
					if(refseq.find(chr)!=refseq.end()){
						cerr<<"Same contig "<<chr<<"in file "<<para.ref<<endl;
					}else{
						refseq.insert(make_pair(chr,seq));
						seqlen.insert(make_pair(chr,seq.length()));
					}
					seq="";
				}
				bin=split(instr);
				string tmpchr=bin[0].substr(1,bin[0].size());
				chr=chr_check(tmpchr);
			}else{
				seq+=instr;
			}
		}
		if(refseq.find(chr)!=refseq.end()){
			cerr<<"Same contig "<<chr<<"in file "<<para.ref<<endl;
		}else{
			refseq.insert(make_pair(chr,seq));
			seqlen.insert(make_pair(chr,seq.length()));
		}
		seq="";
}

void strtoivec(vector<string> &svec, vector<int> &ivec)
{
	for (size_t i=0; i<svec.size();i++)
	{
		ivec.push_back(atoi(svec[i].c_str()));
	}
}

inline string inTstring(int &in)
{
	string s;
	stringstream ss;
	ss<<in;
	ss>>s;
	return s;
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

// A function to spilt string s into vector vec according to char splitchar
void StringSplit(std::string s, char splitchar, std::vector<std::string>& vec) {
	// promise the vec is empty
	if(vec.size()>0)
		vec.clear();
	// chomp the string
	while(s.size() > 0 && s[0] == splitchar)
		s = s.substr(1, s.length() - 1);
	while(s.size() > 0 && s[s.length() - 1] == splitchar)
		s = s.substr(0, s.length() - 1);
	int length = s.length();
	int start=0;
	for(int i=0;i<length;i++)
	{
		if(s[i] == splitchar)
		{
			vec.push_back(s.substr(start,i - start));
			while(s[i + 1] == splitchar)
				i ++;
			start = i+1;
		}
		else if(i == length-1)// attach last
		{
			vec.push_back(s.substr(start,i+1 - start));
		}
	}
}

void usage(){
	cerr << endl;
	cerr << "BGI2294  "<<PROGRAM_COMPILE_DATE<<' '<<PROGRAM_COMPILE_TIME<<endl;
	cerr << endl;
	cerr << "Usage:  " <<"MutDet -i in.bam -r ref.fa -q 20 -Q 20 -v mism.txt -o mut.txt.gz -u" << endl;
	cerr << '\t' << "-i [str] input sorted alignment file, required" << endl;
	cerr << '\t' << "-m [str] format of input file, 'sam' for sam, 'bam' for bam, default[bam]" << endl;
	cerr << '\t' << "-r [str] reference fasta file" << endl;
	cerr << '\t' << "-q [int] minimal mapping quality, default[20]" << endl;
	cerr << '\t' << "-Q [int] minimal base quality, default[20]" << endl;
	cerr << '\t' << "-s [int] base quality shift, 64 for Illumina reads, 33 for Sanger reads, default[33]" << endl;
	cerr << '\t' << "-d [int] minimal effective depth required, default[4]" << endl;
	cerr << '\t' << "-a [int] minimal effective alternative allele supporting reads, defalt[2]" << endl;
	cerr << '\t' << "-l [int] offsets from the starts of the reads, default[5]" << endl;
	cerr << '\t' << "-n [int] offsets from the ends of the reads, default[10]" << endl;
	cerr << '\t' << "-e [float] no more than float would be due to the error rate, default[0.01]" << endl;
	cerr << '\t' << "-v [float] file including the backgroud mismatch rate, reqiured" << endl;
	cerr << '\t' << "-o [str] output file [gz format], required" << endl;
	cerr << '\t' << "-u discard the non-unique mapping reads, default[off]" << endl;
    cerr << '\t' << "-t emit all sites, default[off]" << endl;
	cerr << '\t' << "-h print this usage\n" << endl;
    cerr << "Note: the file of -v can be produced by program MismatchStat, or should include one line format as:\nMismatch rate     : 0.00864522" << endl;
	cerr << endl;
	return;
}

