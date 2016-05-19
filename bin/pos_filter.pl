#!/usr/bin/perl
use strict;
use warnings;
#===============================================================================
#       AUTHOR:  Dongbing Liu, liudongbing2@genomics.org.cn
#      CREATED:  2011/09/08
#     MODIFIED:  2011/10/25 add the --form parameter to sepecify the columns of chr and pos in input file
#===============================================================================
use Getopt::Long;
use File::Basename;
my $usage=<<USAGE;
Description
    This program is used to filter the --input file by using sites in --pos or --loci file, filter means keep(--filt 1) or filter out(--filt 2); of the two kinds of sites file, at least one should be setted! If --pos setted, then the input file should only include one chromosome, otherwise program will go wrong!
Parameter
	--input	[STR] input file to filter, should have chr and pos information
	--form  [STR] num1:num2, column of chr and pos of input file, default[1:2]
	--pos	[STR] a position file with one column: pos, default NULL
	--loci	[STR] a loci file with two columns: chr pos, default NULL
	--filt	[INT] 1|2, 1 represents keep the locis which are in --pos or --loci file, while 2 means filter out, default[1]
	--out	[STR] output file
	--oform	[INT] 1|2, 1 represents output the whole line, while 2 means output only: chr pos, default[1]
	--gz	means output gzip compressed file
	--help	give this information
Exmple 
	perl $0 --input chr.vcf.gz --pos pos.txt --filt 1 --out chr.flt.vcf.gz --oform 2 --gz
	perl $0 --input chr.vcf.gz --loci loci.txt --filt 1 --out chr.flt.vcf.gz --oform 1 --gz --form 4:5
USAGE
my ($input,$form,$posfile,$locifile,$filter,$out,$oform,$gz,$help);
GetOptions(
    "input=s"=>\$input,
    "form=s"=>\$form,
	"pos=s"=>\$posfile,
	"loci=s"=>\$locifile,
	"filter=i"=>\$filter,
    "out=s"=>\$out,
	"oform=i"=>\$oform,
	"gz"=>\$gz,
    "help" => \$help,
);
die "$usage" if(!$input || !$out || $help);
die "$usage" if(!$posfile && !$locifile);
$form||="1:2";
$filter||=1;
$oform||=1;
my ($chr_col,$pos_col);
$form=~/:/ ? ($chr_col,$pos_col)=split /:/,$form : die "--form should have the right as the sample\n";
my $chr_idx=$chr_col-1;
my $pos_idx=$pos_col-1;

print `date`;
my %pos=();
my %loci=();
my $pos_num=0;
if($posfile){
	print "Reading $posfile\n";
	if($posfile=~/\.gz$/){
		open POS,"<:gzip","$posfile" or die $!;
	}else{
		open POS,"$posfile" or die $!;
	}
	while(<POS>){
		chomp;
		if($_=~/^\#/){next;}
		$pos{$_}=1;
        $pos_num++;
	}
	close POS;
}
my $loci_num=0;
if($locifile){
	print "Reading $locifile\n" or die $!;
	if($locifile=~/\.gz$/){
		open LOCI,"gzip -dc $locifile | " or die $!;
	}else{
		open LOCI,"$locifile" or die $!;
	}
	while(<LOCI>){
		chomp;
		if($_=~/^\#/){next;}
		my @loci_line=split /\s+/,$_;
        $loci_line[0]=~s/^chr//;
		$loci{$loci_line[0]}{$loci_line[1]}=1;
        $loci_num++;
	}
	close LOCI;
}
if($input=~/\.gz$/){
	open IN,"gzip -dc $input | " or die $!;
}else{
	open IN,"$input" or die $!;
}
if($gz){
	open OUT,">:gzip","$out" or die $!;
}else{
	open OUT,">$out" or die $!;
}
my $m=0;
my $n=0;
print "Reading $input\n" or die $!;
while(<IN>){
	chomp;
	if($_=~/^\#/){if($oform==1){print OUT "$_\n";}next;}
	$m++;
	my @line=split /\s+/,$_;
    my $chr=$line[$chr_idx];
    my $pos=$line[$pos_idx];
    $chr=~s/^chr//;
	if($posfile){
		if($filter==1){
			next unless (exists $pos{$pos});
		}elsif($filter==2){
			if(exists $pos{$pos}){next;}
		}
	}
	if($locifile){
		if($filter==1){
			next unless (exists $loci{$chr}{$pos});
		}elsif($filter==2){
			if(exists $loci{$chr}{$pos}){next;}
		}
	}
	$n++;
	if($oform==1){print OUT "$_\n";}elsif($oform==2){print OUT "$line[$chr_idx]\t$pos\n";}
}
close IN;
close OUT;
if($posfile){print "Num of sites in pos file : $pos_num\n";}
if($locifile){print "Num of sites in loci file : $loci_num\n";}
print "Total: $m\nLeft: $n\n";
print `date`;
