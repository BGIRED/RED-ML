# RED-ML
RNA editing detection based on machine learning

Description:

    This program is used to identify RNA editing sites based on machine learning.  All of the programs were located in directory of bin when the software package was decompressed. The software package could only run on linux platform currently and the main program is red_ML.pl.

Parameters:

    --rnabam       [STR] the RNA alignment file which is used to calling RNA editing sites.
    --reference    [STR] the genome fasta file(hg19.fa).
    --dbsnp        [STR] SNP database, eg dbsnp138.
    --simpleRepeat [STR] genome simple repeat region annotation file, should be bed format.
    --alu          [STR] genome alu region annotation file, should be bed format.
    --snplist      [STR] a file included konwn SNP sites, the first two column should be chromosome, position and seperated by Tab, it could be used to remove SNP [optional].
    --outdir       [STR] the directory of output.
    --p            [NUM] a number between 0 and 1 [default 0.5];
    --help         [STR] show this information!

Example:

    perl red_ML.update.pl --rnabam in.bam --reference hg19.fa --dbsnp dbsnp138.vcf --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --snplist snp.list --outdir outdir
    Besides, we have offered the example datasets in the software package. You could test the example to guarantee the software package you have downloaded to be completed and it need only about 2 minutes to run the example compleletly.

Requirements:

    Before running this program, several pieces of database are needed.
    dbSNP138, downloaded from:http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
    simpleRepeat, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
    awk '{print $2"\t"$3"\t"$4}' simpleRepeat.txt > simpleRepeat.bed
    bedtools merge -i simpleRepeat.bed > simpleRepeat.merge.bed
    Alu, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
    grep Alu rmsk.txt | awk '{print $6"\t"$7"\t"$8}' > hg19.alu.bed
    We have handled the required databases, and it offered in software package. You could downloaded all the databases from githup, or handle it by yourself in other way. 

Outputs:

    hen the program running completed, three files would be produced in the output directory.
    RNA_editing.sites.txt is the all reliable RNA editing sites file;
    variation.sites.feature.txt is the all varitaion sites file which contain feature values;
    mut.txt.gz is the all varitaion sites file which contain pileup information;

Notice:

    The input bam should be indexed, you could use samtools to create index.
    samtools index in.bam

