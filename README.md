# RED-ML: RNA editing detection based on machine learning

## Description

RED-ML is a software tool to do genome-wide RNA editing dectection (RED) based on RNA-seq data. All source codes and executables are located in the "bin" directory. The tool can be run on a Linux platform and the main program is red_ML.pl.

## Parameters

    --rnabam       [STR] the sorted BAM file obtained from RNA-seq to detect RNA editing sites.
    --reference    [STR] the fasta file containing the reference genome, e.g., hg19.fa.
    --dbsnp        [STR] the SNP database file, e.g., dbSNP138.
    --simpleRepeat [STR] genome-wide simple repeat annotation, should be in BED format.
    --alu          [STR] genome-wide Alu repeat annotation, should be in BED format.
    --snplist      [STR] a tab-delimited file listing known SNPs, with the first two columns being chromosome and position of each SNP [optional].
    --outdir       [STR] the directory of output.
    --p            [NUM] the detection threshold, a number between 0 and 1 [default 0.5];
    --help         [STR] show this help information!

## Example

We have provided a simple example to test the installation of RED-ML. Under the "example" directory, run:
   Xiong Heng, please fill in the command here
It should finish running in ~2 minutes with the following output files (again, please fill in details). Here is another example of using RED-ML:
    perl red_ML.update.pl --rnabam in.bam --reference hg19.fa --dbsnp dbsnp138.vcf --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --snplist snp.list --outdir outdir

## Requirements

RED-ML requires the following data files at the time of public release:
- The reference genome (hg19), downloaded from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes.
- dbSNP138, downloaded from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database.
- simpleRepeat, downloaded from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database, and then do:

    awk '{print $2"\t"$3"\t"$4}' simpleRepeat.txt > simpleRepeat.bed
    bedtools merge -i simpleRepeat.bed > simpleRepeat.merge.bed
    
- Alu, downloaded from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database, and do:

    grep Alu rmsk.txt | awk '{print $6"\t"$7"\t"$8}' > hg19.alu.bed

We have also provided the simpleRepeat and Alu files under the "database" directory for the user's convenience.

## Optional

### SNP calling

If you have matching DNA-seq data or aligned DNA BAM files, we strongly recommend to take advantage of them. You could call SNPs by GATK (haplotypecaller) or SOAPsnp and modify the format of the resulting file (such as vcf) to fit the format required by --snplist.
    
### Alignment

Although RED-ML can accept BAM files produced by different alignment tools, the current version has been optimized for BWA and TopHat2 during the construction of our ML model, and we find that the choice of alignment tools and the associated parameters could have a large impact on RED. To help users with proper alignment strategies, we recommend the following steps:

1. When reads are aligned by BWA (preferred), one should first build a new reference which combines reference genome (hg19) and exonic sequences surrounding all known splice junctions, and the detail method is the same as in Ramaswami et al. (Nature Methods 2012) and Wang et al. (GigaScience 2016). SAMtools can be used to sort the alignment file and remove the PCR duplicate reads.

2. When TopHat2 is chosen, the cleaned reads can be mapped to the reference genome (hg19) directly with default parameters. Picard should be used to sort the alignment and to remove duplicate reads induced by PCR, and base quality score recalibration can be carried out by GATK.
    
## Outputs

When the program finishes running, three files will be created in the output directory. RNA_editing.sites.txt lists all detected RNA editing sites that pass the detection threshold p; variation.sites.feature.txt lists all variant sites with associated feature values; mut.txt.gz contains all variant sites with pileup information.

## Notice

The input bam should be sorted (indexed), you could use samtools to create index.

    samtools index in.bam

