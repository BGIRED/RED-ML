# RED-ML
RNA editing detection based on machine learning

Description:

    RED-ML is a software tool to do genome-wide RNA editing dectection (RED) based on RNA-seq data.
    All source codes and executables are located in the "bin" directory.
    The tool can be run on a Linux platform and the main program is red_ML.pl.

Parameters:

    --rnabam       [STR] the RNA alignment file which is used to calling RNA editing sites.
    --reference    [STR] the genome fasta file (hg19.fa).
    --dbsnp        [STR] SNP database, eg dbSNP138.
    --simpleRepeat [STR] genome simple repeat region annotation file, should be BED format.
    --alu          [STR] genome alu region annotation file, should be bed format.
    --snplist      [STR] a file included konwn SNP sites, the first two column should be chromosome, position and seperated by Tab, it could be used to remove SNP [optional].
    --outdir       [STR] the directory of output.
    --p            [NUM] a number range from 0 to 1 [default 0.5];
    --help         [STR] show this information!

Example:

    perl red_ML.update.pl --rnabam in.bam --reference hg19.fa --dbsnp dbsnp138.vcf --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --snplist snp.list --outdir outdir
    We have also provided a simple example to test the installation under the "example" directory.
    The command to run the example is: (Xiong Heng, please fill it in).
    It should complete in ~2 minutes.

Requirements:

    RED-ML requires the following data files at the time of public release:
    The reference genome (hg19), downloaded from : http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes.
    dbSNP138, downloaded from:http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
    simpleRepeat, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database, and then do:
    	awk '{print $2"\t"$3"\t"$4}' simpleRepeat.txt > simpleRepeat.bed
    	bedtools merge -i simpleRepeat.bed > simpleRepeat.merge.bed
    Alu, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database, and do:
    	 grep Alu rmsk.txt | awk '{print $6"\t"$7"\t"$8}' > hg19.alu.bed
    We have also provided the simpleRepeat and Alu files under the "database" directory.

optional:

    SNP calling:
    If you have DNA datasets, such as DNA alignment, we strongly recommend to take advantage of it. You could call SNPs by GATK (haplotypecaller) or SOAPsnp and moidfy the format of the resulting file (such as vcf) to fit the format required by --snplist.
    
    Alignment:
    Although RED-ML can accept BAM files produced by different alignment tools, the current version has been optimized for BWA and TopHat2 during the construction of the ML model, and we find that the choice of alignment tools and the associated parameters could have a large impact on RED. To help users with proper alignment strategies, we recommend the following steps:
    (1) When reads are aligned by BWA (preferred), one should first build a new reference which combines reference genome (hg19) and exonic sequences surrounding all known splice junctions, and the detail method is the same as in Ramaswami et al. (Nature Methods 2012) and Wang et al. (GigaScience 2016). SAMtools can be used to sort the alignment file and remove the PCR duplicate reads.
    (2) When TopHat2 is chosen, the cleaned reads can be mapped to the reference genome (hg19) directly with default parameters. Picard should be used to sort the alignment and to remove duplicate reads induced by PCR, and base quality score recalibration can be carried out by GATK.
   
    
Outputs:

    When the program finishes running, three files will be created in the output directory.
    RNA_editing.sites.txt lists all detected RNA editing sites that pass the detection threshold p;
    variation.sites.feature.txt lists all variant sites with associated feature values;
    mut.txt.gz contains all variant sites with pileup information.

Notice:

    The input bam should be sorted (indexed), you could use samtools to create index.
    samtools index in.bam

