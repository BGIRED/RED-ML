# RED-ML
RNA editing detection based on machine learning

Description:

    RED-ML is a software tool to do genome-wide RNA editing dectection (RED) based on RNA-seq data. All source codes and executables are located in the "bin" directory.  The tool can be run on a Linux platform and the main program is red_ML.pl.

Parameters:

    --rnabam       [STR] the RNA alignment file which is used to calling RNA editing sites.
    --reference    [STR] the genome fasta file(hg19.fa).
    --dbsnp        [STR] SNP database, eg dbsnp138.
    --simpleRepeat [STR] genome simple repeat region annotation file, should be bed format.
    --alu          [STR] genome alu region annotation file, should be bed format.
    --snplist      [STR] a file included konwn SNP sites, the first two column should be chromosome, position and seperated by Tab, it could be used to remove SNP [optional].
    --outdir       [STR] the directory of output.
    --p            [NUM] a number range from 0 to 1 [default 0.5];
    --help         [STR] show this information!

Example:

    perl red_ML.update.pl --rnabam in.bam --reference hg19.fa --dbsnp dbsnp138.vcf --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --snplist snp.list --outdir outdir
    Besides, we have offered the example datasets in the software package. You could test the example to guarantee the software package you have downloaded to be completed and it need only about 2 minutes to run the example compleletly.

Requirements:

    Before running this program, several pieces of database are needed.
    reference (hg19), downloaded from : http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes.
    dbSNP138, downloaded from:http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
    simpleRepeat, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
    awk '{print $2"\t"$3"\t"$4}' simpleRepeat.txt > simpleRepeat.bed
    bedtools merge -i simpleRepeat.bed > simpleRepeat.merge.bed
    Alu, downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
    grep Alu rmsk.txt | awk '{print $6"\t"$7"\t"$8}' > hg19.alu.bed
    We have handled the required databases, and it offered in software package. You could downloaded all the databases from githup, or handle it by yourself in other way. 

optional:

    SNP calling:
    If you have DNA datasets, such as DNA alignment, we strongly recommend to take advantage of it. You could call SNP based on this datasets by GATK (haplotypecaller) or SOAPsnp and moidfy the format of result file (such as vcf) to parameter --snplist demanded format.
    
    Alignment:
    Although RED-ML can accept BAM files produced by different alignment tools, the current version has been specifically optimized for BWA and TopHat2 due to the construction of model, and we find that the choice of alignment tools and the associated parameters could have a large impact on RED. To help users with proper alignment strategies, we have detailed some recommendations:
    When reads are aligned by BWA (recommended), one should first build a new reference which combines reference genome (hg19) and exonic sequences surrounding all known splice junctions, and the detail method is same as Ramaswami et al. and Wang et al. SAMtools was used to sort the alignment file and to remove the PCR duplicate reads.
    When TopHat2 a fast splice junction mapper for RNA-seq reads is selected, the cleaned reads were mapped to reference genome (hg19) directly with default parameters. Picard was used to sort the alignment and to remove duplicate reads induced by PCR, and base quality score recalibration was carried out by GATK.
   
    
Outputs:

    When the program finishes running, three files will be created in the output directory.
    RNA_editing.sites.txt lists all detected RNA editing sites that pass the detection threshold p;
    variation.sites.feature.txt lists all variant sites with associated feature values;
    mut.txt.gz contains all variant sites with pileup information.

Notice:

    The input bam should be sorted (indexed), you could use samtools to create index.
    samtools index in.bam

