#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;

my $usage=<<USAGE;

     Description:

        RED-ML is a software tool to do genome-wide RNA editing dectection (RED) based on RNA-seq data. All source codes and executables are located in the "bin" directory. The tool can be run on a Linux platform and the main program is red_ML.pl.
    
    Parameters:
        
        --rnabam       [STR] the sorted BAM file obtained from RNA-seq to detect RNA editing sites.
        --reference    [STR] the fasta file containing the reference genome, e.g., hg19.fa.
        --dbsnp        [STR] the SNP database file, e.g., dbSNP138.
        --simpleRepeat [STR] genome-wide simple repeat annotation, should be in BED format.
        --alu          [STR] genome-wide Alu repeat annotation, should be in BED format.
        --snplist      [STR] a tab-delimited file listing known SNPs, with the first two columns being chromosome and position of each SNP [optional].
        --outdir       [STR] the directory of output.
        --p            [NUM] the detection threshold, a number between 0 and 1 [default 0.5];
        --help         [STR] show this help information!
    
    Example:
    
        perl $0 --rnabam in.bam --reference hg19.fa --dbsnp dbsnp138.vcf --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --snplist snp.list --outdir outdir
     
    Requirements:

        Before running this program, several pieces of database are needed.
        The reference genome (hg19), downloaded from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes.
        dbSNP138, downloaded from:http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
        simpleRepeat, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database, and then do:
        awk '{print \$2"\\t"\$3"\\t"\$4}' simpleRepeat.txt > simpleRepeat.bed
        bedtools merge -i simpleRepeat.bed > simpleRepeat.merge.bed
        Alu, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database, and then do:
        grep Alu rmsk.txt | awk '{print \$6"\\t"\$7"\\t"\$8}' > hg19.alu.bed
        We have also provided the simpleRepeat and Alu files under the "database" directory for the user's convenience.

    Outputs:
 
        When the program finishes running, three files will be created in the output directory. 
        RNA_editing.sites.txt lists all detected RNA editing sites that pass the detection threshold p;
        variation.sites.feature.txt lists all variant sites with associated feature values; 
        mut.txt.gz contains all variant sites with pileup information.
    
    Notice:
    
        The input bam should be sorted (indexed), you could use samtools to create index.
        samtools index in.bam

USAGE
         
my ($rnabam,$reference,$dbsnp,$simpleRepeat,$alu,$snplist,$outdir,$help,$p);
GetOptions(
    "rnabam=s"=>\$rnabam,
    "reference=s"=>\$reference,
    "dbsnp=s"=>\$dbsnp,
    "simpleRepeat=s"=>\$simpleRepeat,
    "alu=s"=>\$alu,
    "snplist=s"=>\$snplist,
    "outdir=s"=>\$outdir,
    "p=f"=>\$p,
    "help" => \$help,
);
die "$usage" if(!$rnabam || $help || !$outdir);

print "Begin time: ".`date`;

unless (-d $outdir){`mkdir -p $outdir`;}

$p ||= 0.5;
#Coefficients of features;
my ($c1,$c2,$c3,$c4,$c5,$c6,$c7,$c8,$c9,$c10,$c11,$c12,$c13,$c14,$c15,$c16,$c17,$c18,$c19,$c20,$c21,$c22,$c23,$c24,$c25,$c26,$c27,$c28,$c29) = (0.20814042,28.15978658,29.46200811,0.10715825,6.91274889,-0.50766250,-0.29390003,1.40012066,-3.57552459,0.01493510,0.52292404,-0.06624730,-1.14912825,-0.69362344,1.57646089,-0.00591220,-1.44597477,0.27586286,-0.11703504,-0.08778774,-2.05264373,-0.01213383,-0.05561199,7.48881702,1.20189923,0.37863824,-6.47576116,-5.42376070,-43.72477118);

#Calculate mismatch ratio;
print "Calculating mismatch ratio...\n";
`$Bin/MismatchStat -i $rnabam -o $outdir/stat.txt -u -q 20`;

#Calling varitation sites;
print "Variation sites detecting...\n";
`$Bin/MutDetML -i $rnabam -r $reference -v $outdir/stat.txt -u -o $outdir/mut.txt.gz`;

#Filtering SNP using known SNP sites;
if (defined $snplist){
    print "Filtering SNP...\n";
    `perl $Bin/pos_filter.pl --input $outdir/mut.txt.gz --loci $snplist --filt 2 --out $outdir/mut.snp.txt.gz --gz`;
    `mv $outdir/mut.snp.txt.gz $outdir/mut.txt.gz`;
}

#SNP annotation using dbsnp;
print "SNP annotating...\n";
`perl $Bin/pos_filter.pl --input $dbsnp --loci $outdir/mut.txt.gz --out $outdir/snp.txt`;
my %hashd;
open SNP, "$outdir/snp.txt" or die $!;
while (<SNP>){
    chomp;
    if ($_ =~ /^#/) {next;}
    my @array = split /\s+/,$_;
    my $snpsite = $array[0]."_".$array[1];
    $hashd{$snpsite} = 1;
}
close SNP;

#Store genome position information;
print "Storing genome position information...\n";
open REF, $reference or die $!;
$/ = ">";<REF>;
my %hashr;
while (<REF>){
    chomp;
    my @array = split /\n/,$_;
    my $chr = (split /\s+/,$array[0])[0];
    if ($chr !~ /chr/){$chr = "chr".$chr;}
        shift(@array);
        my $seq = join("",@array);
        $hashr{$chr} = $seq;
}
close REF;
$/ = "\n";

#Simple repeat annotation using simplerepeat database;
print "Simple repeat annotating...\n";
`perl $Bin/pos_filterByRegion.pl $simpleRepeat $outdir/mut.txt.gz > $outdir/out.rep`;
open SR, "$outdir/out.rep" or die $!;
my %hashs;
while (<SR>){
    chomp;
    my @array = split /\s+/,$_;
    my $pos = $array[0]."_".$array[1];
    $hashs{$pos} = 1;
}
close SR;
`rm -f $outdir/out.rep`;

#Alu region annotation using Alu database;
print "Alu annotating...\n";
`perl $Bin/pos_filterByRegion.pl $alu $outdir/mut.txt.gz > $outdir/alu.txt`;
open ALU, "$outdir/alu.txt" or die $!;
my %hasha;
while (<ALU>){
    chomp;
    my @array = split /\s+/,$_;
    my $pos = $array[0]."_".$array[1];
    $hasha{$pos} = 1;
}
close ALU;
`rm -f $outdir/alu.txt`;

#Feature statistics;
print "Feature values calculating...\n";
open IN, "gzip -dc $outdir/mut.txt.gz |" or die $!;
my $first = <IN>;
chomp  $first;
open OUT1, ">$outdir/variation.sites.feature.txt" or die $!;
open OUT2, ">$outdir/RNA_editing.sites.txt" or die $!;
print OUT1 "#Pos\tVariation\tFrequency\tmfre\tQuality\tBino\tRef_end\tRef_mid\tAlt_end\tend_c\tAlt_mid\tend_p\tEndratio\tRef_minus\tRef_plus\tAlt_minus\tAlt_plus\tsb_c\tStrandbias_p\tStrandbiasratio\tHomolen\tSimplerepeat\tLeftvalue\tRightvalue\tAG\tAlu\tmotif\tDBSNP\thighfre"."\n";
print OUT2 "#P_edit: The probability of being a RNA editing site predicted by RED-ML\n";
print OUT2 "#Chromosome\tPosition\tRead_depth\tReference\tReference_support_reads\tAlternative\tAlternative_support_reads\tP_edit\n";
while (<IN>){
    chomp;
    my ($var_reads,$fre,$mfre,$bqua,$bino,$ref_end,$ref_mid,$alt_end,$end_c,$alt_mid,$end_p,$endratio,$ref_minus,$ref_plus,$alt_minus,$alt_plus,$sb_c,$strandbias_p,$strandbiasratio,$homolen,$sr,$leftvalue,$rightvalue,$AG,$alu,$motif,$snpvalue,$alt_end_f,$sb_f,$lowfre,$highfre);
    if ($_ =~ /^#/) {next;}
    my @a = split /\s+/,$_;
    my $pos = $a[0]."_".$a[1];
    $var_reads = log($a[13])/log(10);
    $fre = sprintf "%0.3f",($a[13]/$a[2]);
    if ($fre>=0.95){$highfre=1;}
    else {$highfre=0;}
    my $temp = 1 - $fre;
    $mfre = $temp > 0.9 ? 0.9 : $temp;
    $lowfre = $fre < 0.1 ? 1 : 0;
    $bqua = $a[20];
    $bino = $a[13] > $a[21] ? 1 : 0;
    if ($a[6] > 0) {$ref_end =  log($a[6])/log(10);}
    else {$ref_end = -1;}
    if ($a[5] > 0) {$ref_mid =  log($a[5])/log(10);}
    else {$ref_mid = -1;}
    if ($a[15] > 0) {$alt_end =  log($a[15])/log(10);}
    else {$alt_end = -1;}
    $end_c = $a[15]/$a[13];
    $alt_end_f = $end_c > 0.9 ? 1 : 0;
    if ($a[14] > 0) {$alt_mid =  log($a[14])/log(10);}
    else {$alt_mid = -1;}
    $end_p =  $a[22];
    if ($a[6] ==0 || $a[5] ==0 || $a[14] ==0) {$endratio = 10}
    else {
        $endratio = ($a[15]/$a[14])/($a[6]/$a[5]);
        if ($endratio >= 10) {$endratio = 10;}
    }
    if ($a[8] > 0) {$ref_minus =  log($a[8])/log(10);}
    else {$ref_minus = -1;}
    if ($a[7] > 0) {$ref_plus =  log($a[7])/log(10);}
    else {$ref_plus = -1;}
    if ($a[17] > 0) {$alt_minus =  log($a[17])/log(10);}
    else {$alt_minus = -1;}
    if ($a[16] > 0) {$alt_plus =  log($a[16])/log(10);}
    else {$alt_plus = -1;}
    $sb_c = $a[17]/$a[13];
    if ($sb_c > 0.9 || $sb_c < 0.1) {$sb_f = 1;}
    else {$sb_f = 0;}
    $strandbias_p =  $a[23];
    if ($a[7] ==0 || $a[8] ==0 || $a[16] ==0) {$strandbiasratio = 10;}
    else {
        $strandbiasratio = ($a[17]/$a[16])/($a[8]/$a[7]);
        if ($strandbiasratio >=10) {$strandbiasratio =10;}
    }
    my $refbase = substr($hashr{$a[0]},$a[1]-1,1);
    my $lef = $a[1];
    my $righ = $a[1];
    while (1) {
        $lef--;
        my $temp = substr($hashr{$a[0]},$lef-1,1);
        if ($temp ne $refbase) {last;}
    }
    while (1){
        $righ++;
        my $temp = substr($hashr{$a[0]},$righ-1,1);
        if ($temp ne $refbase) {last;}
    }
    $homolen = ($righ-1)-($lef+1) +1;
    $sr = (exists $hashs{$pos}) ? 1 : 0;
    my $leftbase = uc(substr($hashr{$a[0]},$a[1]-1-1,1)); 
    my $rightbase = uc(substr($hashr{$a[0]},$a[1]-1+1,1));
    my $midbase = uc(substr($hashr{$a[0]},$a[1]-1,1));    
    if ($midbase eq "N"){next;}
    my $left = $leftbase.$midbase;                            
    my $right = $midbase.$rightbase;                          
    if ($left eq "AA" ) {$leftvalue = 1;}                     
    if ($left eq "AT" ) {$leftvalue = 2;}                     
    if ($left eq "AC" ) {$leftvalue = 3;}                     
    if ($left eq "AG" ) {$leftvalue = 4;}                     
    if ($left eq "TA" ) {$leftvalue = 5;}                     
    if ($left eq "TT" ) {$leftvalue = 6;}                     
    if ($left eq "TC" ) {$leftvalue = 7;}                     
    if ($left eq "TG" ) {$leftvalue = 8;}                     
    if ($left eq "CA" ) {$leftvalue = 9;}                     
    if ($left eq "CT" ) {$leftvalue = 10;}                    
    if ($left eq "CC" ) {$leftvalue = 11;}                    
    if ($left eq "CG" ) {$leftvalue = 12;}                    
    if ($left eq "GA" ) {$leftvalue = 13;}                    
    if ($left eq "GT" ) {$leftvalue = 14;}                    
    if ($left eq "GC" ) {$leftvalue = 15;}                    
    if ($left eq "GG" ) {$leftvalue = 16;}
    if ($leftbase eq "N") {$leftvalue = 17;}
                                                              
    if ($right eq "AA" ) {$rightvalue = 1;}                   
    if ($right eq "AT" ) {$rightvalue = 2;}                   
    if ($right eq "AC" ) {$rightvalue = 3;}                   
    if ($right eq "AG" ) {$rightvalue = 4;}                   
    if ($right eq "TA" ) {$rightvalue = 5;}                   
    if ($right eq "TT" ) {$rightvalue = 6;}                   
    if ($right eq "TC" ) {$rightvalue = 7;}                   
    if ($right eq "TG" ) {$rightvalue = 8;}                   
    if ($right eq "CA" ) {$rightvalue = 9;}                   
    if ($right eq "CT" ) {$rightvalue = 10;}                  
    if ($right eq "CC" ) {$rightvalue = 11;}                  
    if ($right eq "CG" ) {$rightvalue = 12;}                  
    if ($right eq "GA" ) {$rightvalue = 13;}                  
    if ($right eq "GT" ) {$rightvalue = 14;}                  
    if ($right eq "GC" ) {$rightvalue = 15;}                  
    if ($right eq "GG" ) {$rightvalue = 16;}
    if ($rightbase eq "N" or $rightbase eq "") {$rightvalue = 17;}
    
    if ($a[3] eq "A" and $a[12] eq "G") {$AG = 1;}
    elsif ($a[3] eq "T" and $a[12] eq "C"){$AG = 1;}
    else {$AG = 0;}
    $alu = (exists $hasha{$pos}) ? 1 : 0;
    if ($midbase eq "A" and $rightbase eq "G" and $leftbase ne "G") {$motif = 1;}
    elsif ($midbase eq "T" and $rightbase ne "C" and $leftbase eq "C") {$motif = 1;}
    else {$motif = 0;}
    $snpvalue = (exists $hashd{$pos}) ? 1 : 0;
    print OUT1 $pos."\t".$var_reads."\t".$fre."\t".$mfre."\t".$bqua."\t".$bino."\t".$ref_end."\t".$ref_mid."\t".$alt_end."\t".$end_c."\t".$alt_mid."\t".$end_p."\t".$endratio."\t".$ref_minus."\t".$ref_plus."\t".$alt_minus."\t".$alt_plus."\t".$sb_c."\t".$strandbias_p."\t".$strandbiasratio."\t".$homolen."\t".$sr."\t".$leftvalue."\t".$rightvalue."\t".$AG."\t".$alu."\t".$motif."\t".$snpvalue."\t".$highfre."\n";
    my $Z = -(($c1)*$var_reads+($c2)*$fre+($c3)*$mfre+($c4)*$bqua+($c5)*$bino+($c6)*$ref_end+($c7)*$ref_mid+($c8)*$alt_end+($c9)*$end_c+($c10)*$alt_mid+($c11)*$end_p+($c12)*$endratio+($c13)*$ref_minus+($c14)*$ref_plus+($c15)*$alt_minus+($c16)*$alt_plus+($c17)*$sb_c+($c18)*$strandbias_p+($c19)*$strandbiasratio+($c20)*$homolen+($c21)*$sr+($c22)*$leftvalue+($c23)*$rightvalue+($c24)*$AG+($c25)*$alu+($c26)*$motif+($c27)*$snpvalue+($c28)*$highfre+($c29));
    my $result =  1/(1+2.718281828459**$Z);
    if ($result >= $p){
        print OUT2 $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]."\t".$a[4]."\t".$a[12]."\t".$a[13]."\t".$result."\n";
    }
}
close IN;
close OUT1;
close OUT2;

print "End time: ".`date`;
