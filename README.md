# TotiEVEs

This document is a walkthrough of the methods and code used to analyze Endogenous Toti-Like Viral Elements (ETLVEs) in three planthopper genomes, including *Sogatella furcifera*, *Laodelphax striatellus*, and *Nilaparvata lugens*. 

## 1 - Distribution of ETLVEs in planthopper individuals

### 1.1 - Mapping

To understand the distribution of ETLVEs, we mapped Illumina reads to planthopper genomes using BWA MEM (version 0.7.17). Note that for your analysis, please replace $sampleID to the actual sample names. Genomes of three planthopper species are downloaded from NCBI (GCA_014465815.1 for *Laodelphax striatellus*, GCA_014356515.1 for *Sogatella furcifera*, and GCA_014356525.1 for *Nilaparvata lugens*). 

        bwa mem -M -t 40 -R '@RG\tID:'"$sampleID"'\tSM:'"$sampleID"'\tPL:Illumina' -o $sampleID.sam SfHau_GCA_014356515.1_genomic.fa ./links/$sampleID_1_clean.fq.gz ./links/$sampleID_2_clean.fq.gz
