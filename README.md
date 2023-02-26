# TotiEVEs

This document is a walkthrough of the methods and code used to analyze Endogenous Toti-Like Viral Elements (ETLVEs) in three planthopper genomes, including *Sogatella furcifera*, *Laodelphax striatellus*, and *Nilaparvata lugens*. 

## 1 - Distribution of ETLVEs in planthopper individuals

### 1.1 - Mapping

To understand the distribution of ETLVEs, we mapped Illumina reads to planthopper genomes using BWA MEM (version 0.7.17). Note that for your analysis, please replace $sampleID to the actual sample names. Genomes of three planthopper species are downloaded from NCBI (GCA_014465815.1 for *Laodelphax striatellus*, GCA_014356515.1 for *Sogatella furcifera*, and GCA_014356525.1 for *Nilaparvata lugens*). 

        # BWA MEM
        ~/bin/bwa mem -M -t 40 -R '@RG\tID:'"$sampleID"'\tSM:'"$sampleID"'\tPL:Illumina' -o $sampleID.sam SfHau_GCA_014356515.1_genomic.fa ./links/$sampleID_1_clean.fq.gz ./links/$sampleID_2_clean.fq.gz
        
        # Convert SAM files to BAM files
        ~/bin/samtools sort --threads 40 $sampleID.bam > $sampleID.sorted.bam
        ~/bin/samtools index -@ 40 $sampleID.sorted.bam
        
        # MarkDuplications 
        ~/bin/gatk-4.2.6.1/gatk MarkDuplicates -I $sampleID.sorted.bam -O $sampleID_markdup.bam -M $sampleID_markdup_metrics.txt
        
        # Index BAM files
        ~/bin/samtools sort --threads 40 $sampleID_markdup.bam > $sampleID_markdup.sorted.bam
        ~/bin/samtools index -@ 40 $sampleID_markdup.sorted.bam
        
### 1.2 - SNP calling using GATK haplotypeCaller
        ~/bin/gatk-4.2.6.1/gatk HaplotypeCaller -R NlHau_GCA_014356525.1_genomic.fa -I $sampleID_markdup.sorted.bam -O $sampleID_markdup.EVE.gvcf.gz --native-pair-hmm-threads 40 --intervals EVE_fast_slow.sort.list -ERC GVCF --disable-read-filter MappingQualityReadFilter
