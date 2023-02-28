# Domesticated non-retroviral integrated RNA virus sequences serving host cellular innovations in arthropods

This document is a walkthrough of the methods and code used to analyze Endogenous Toti-Like Viral Elements (ETLVEs) in three planthopper genomes, including *Sogatella furcifera* (GCA_014356515.1), *Laodelphax striatellus* (GCA_014465815.1), and *Nilaparvata lugens* (GCA_014356525.1). 

## 1 - Distribution of ETLVEs in planthopper individuals

To understand the heterozygosity of EVEs, fast-evolving and slow-evolving genes in planthoppers, we first identified fast-evolving and slow-evolving genes gene and related regions in the genome. 

Genomes of three planthopper species were downloaded from NCBI (GCA_014465815.1 for *Laodelphax striatellus*, GCA_014356515.1 for *Sogatella furcifera*, and GCA_014356525.1 for *Nilaparvata lugens*). 

Genome annotations of the three species were downloaded from 

### 1.1 - EVE loci in planthopper genomes

        # BLAST. BLAST inputs for Sogatella furcifera (EVEs_sf.fa), Laodelphax striatellus (EVEs_ls.fa), and Nilaparvata lugens (EVEs_nl.fa) are available under ./step1.1_BLAST/. Genome sequences were downloaded from NCBI. 
        # Sogatella furcifera
        blastn -db sf_GCA_014356515.1_genomic.fa -query EVEs_sf.fa -out EVEs_sf.blast -evalue 1e-10 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        
        # Laodelphax striatellus
        blastn -db ls_GCA_014465815.1_genomic.fa -query EVEs_ls.fa -out EVEs_ls.blast -evalue 1e-10 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        
        # Nilaparvata lugens
        blastn -db nl_GCA_014356525.1_genomic.fa -query EVEs_nl.fa -out EVEs_nl.blast -evalue 1e-10 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        
        # Only matches with >99% global identity are kept for downstream analyses
        # Filtered BLAST hits can be found here: ./step1.1_BLAST/*.blast.filtered

### 1.2 - Orthologs assignment

        cd ./step1.2_orthology/
        # Find the longest isoform. Protein sequences (Nilaparvata_lugens.anno.pep.fa, Sogatella_furcifera.anno.pep.fa, and Laodelphax_striatellus.anno.pep.fa) were downloaded from http://v2.insect-genome.com/Organism/572 for Nilaparvata lugens, http://v2.insect-genome.com/Organism/709 for Sogatella furcifera, and http://v2.insect-genome.com/Organism/477 for Laodelphax striatellus. 
        perl find_longest_protein_Braker_v20210403b.pl Nilaparvata_lugens.anno.pep.fa Nilaparvata_lugens.anno.pep.fa.longest
        perl find_longest_protein_Braker_v20210403b.pl Sogatella_furcifera.anno.pep.fa Sogatella_furcifera.anno.pep.fa.longest
        perl find_longest_protein_Braker_v20210403b.pl Laodelphax_striatellus.anno.pep.fa Laodelphax_striatellus.anno.pep.fa.longest
        
        # Remove the proteins with stop codons
        less -S Nilaparvata_lugens.anno.pep.fa.longest | perl -e '$/=">"; <>; while(<>){ chomp; @line = split/\n/; $name = shift @line; $seq = join "", @line; if($seq=~/X/){ }elsif(/\*/){ }else{ print ">$name\n$seq\n"; } } ' > Nilaparvata_lugens.anno.pep.fa.noStopCodon
        less -S Sogatella_furcifera.anno.pep.fa.longest | perl -e '$/=">"; <>; while(<>){ chomp; @line = split/\n/; $name = shift @line; $seq = join "", @line; if($seq=~/X/){ }elsif(/\*/){ }else{ print ">$name\n$seq\n"; } } ' > Sogatella_furcifera.anno.pep.fa.longest.noStopCodon
        less -S Laodelphax_striatellus.anno.pep.fa.longest | perl -e '$/=">"; <>; while(<>){ chomp; @line = split/\n/; $name = shift @line; $seq = join "", @line; if($seq=~/X/){ }elsif(/\*/){ }else{ print ">$name\n$seq\n"; } } ' > Laodelphax_striatellus.anno.pep.fa.noStopCodon
        
        # Copy proteins to a new folder for OrthoFinder
        cp Nilaparvata_lugens.anno.pep.fa run1/Nilu.faa
        cp Sogatella_furcifera.anno.pep.fa.longest run1/Sofu.faa
        cp Laodelphax_striatellus.anno.pep.fa run1/Last.faa
        
        # Run OrthoFinder version 2.5.4
        python ~/bin/orthoFinder_v2.5.4/OrthoFinder_source/orthofinder.py -t 40 -a 40 -f run1
        
        # Find 1:1:1 orthologs among three planthopper species
        cp ../run2_noStopCodon/OrthoFinder/Results_Aug09/Orthogroups/Orthogroups.GeneCount.tsv ./
        cp -r ../run2_noStopCodon/OrthoFinder/Results_Aug09/Orthogroup_Sequences/ ./
        
        perl ./step1.2_orthology/filter_orthologs.pl Orthogroups.GeneCount.tsv 1 Orthogroups.GeneCount.tsv.1.0
        
        # Batch mafft-linsi (version 7.505) and Gblocks (version 0.91b)
        perl ./step1.2_orthology/concatenate.pl Orthogroups.GeneCount.tsv.1.0 Orthogroup_Sequences/ Orthogroups.GeneCount.tsv.1.0.fasta


### 1.3 - Search for top 5% fast-evolving and top 5% slow-evolving genes in 1:1:1 single-copy orthologs among three planthopper species

        # Calculate p-distance and gap ratio in batch. Sequences with gap ratio >=10% were removed from downstream analyses. 
        perl ./step1.3_evolutionary_rate/calculate_distance_batch.p_distance.pl Orthogroups.GeneCount.tsv.1.0 ./Ortholog_Alignment/
        
        # Combine all the distance outputs
        grep "Average" ./Ortholog_Alignment/*distance | perl -ne 'if(/\/(\w+)\./){ @a=split; print "$1\t$a[-1]\n"; }' > dist_noGBlocks_gapFiltered.summary
        
        # Assign slow-evolving and fast-evolving genes; we have 1,462 genes, 5% would be 1,462 x 0.05 = 73
        # Therefore, the top 5% and bottom 5% genes were assigned as fast-evolving and slow-evolving genes
        sort -k2,2nr dist_noGBlocks_gapFiltered.summary | head -73 > dist.summary.fastEvolving
        sort -k2,2nr dist_noGBlocks_gapFiltered.summary | tail -73 > dist.summary.slowEvolving
        
        # Find fast-evolving and slow-evolving genes in Sogatella furcifera, Laodelphax striatellus, and Nilaparvata lugens 
        less -S dist.summary.slowEvolving | perl -ne '@a=split; print "$a[0]\n"; `cp ./Ortholog_Alignment/$a[0].GBlocks.fas ../../step6_LoP/slowGenes/`; '

### 1.4 - Fast-evolving and slow-evolving gene loci in planthopper genomes
Fast- and slow-evolving protein sequences were first extracted from orthologous groups for each planthopper species. TBLASTN was then used to align proteins to the corresponding genome. Loci of the genes were recorded for level of polymorphism (LoP) calculation. 

        # Copy genes to fast- or slow-evolving gene folder
        less -S dist.summary.fastEvolving | perl -ne '@a=split; print "$a[0]\n"; `cp ./Orthogroup_Sequences/$a[0].fa ./fastGenes/`; '
        less -S dist.summary.slowEvolving | perl -ne '@a=split; print "$a[0]\n"; `cp ./Orthogroup_Sequences/$a[0].fa ./slowGenes/`; '
        
        # Save protein sequences
        # Sogatella furcifera
        cat ./fastGenes/*.fa | grep -A 1 Sfur | grep -v "^\-\-" > sf.fastGenes.aa.fasta
        cat ./slowGenes/*.fa | grep -A 1 Sfur | grep -v "^\-\-" > sf.slowGenes.aa.fasta
        
        # Laodelphax striatellus
        cat ./fastGenes/*.fa | grep -A 1 Lstr | grep -v "^\-\-" > ls.fastGenes.aa.fasta
        cat ./slowGenes/*.fa | grep -A 1 Lstr | grep -v "^\-\-" > ls.slowGenes.aa.fasta
        
        tblastn -db ls_GCA_014465815.1_genomic.fa -query ls.slowGenes.aa.fasta -out ls.slowGenes.aa.fasta.blast -evalue 1e-5 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 8
        tblastn -db ls_GCA_014465815.1_genomic.fa -query ls.fastGenes.aa.fasta -out ls.fastGenes.aa.fasta.blast -evalue 1e-5 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 8
        
        less -S ls.fastGenes.aa.fasta.blast | awk '($5 >= 80){print }' > ls.fastGenes.aa.fasta.blast.filter
        less -S ls.slowGenes.aa.fasta.blast | awk '($5 >= 80){print }' > ls.slowGenes.aa.fasta.blast.filter

        # Nilaparvata lugens
        cat ./fastGenes/*.fa | grep -A 1 Nlug | grep -v "^\-\-" > nl.fastGenes.aa.fasta
        cat ./slowGenes/*.fa | grep -A 1 Nlug | grep -v "^\-\-" > nl.slowGenes.aa.fasta
        
        tblastn -db nl_GCA_014356525.1_genomic.fa -query nl.fastGenes.aa.fasta -out nl.fastGenes.aa.fasta.blast -evalue 1e-5 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 8
        tblastn -db nl_GCA_014356525.1_genomic.fa -query nl.slowGenes.aa.fasta -out nl.slowGenes.aa.fasta.blast -evalue 1e-5 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 8
        
        less -S nl.fastGenes.aa.fasta.blast | awk '($5 >= 80){print;}' > nl.fastGenes.aa.fasta.blast.filter
        less -S nl.slowGenes.aa.fasta.blast | awk '($5 >= 80){print;}' > nl.slowGenes.aa.fasta.blast.filter

        # Sogatella furcifera
        cat ./fastGenes/*.fa | grep -A 1 Sfur | grep -v "^\-\-" > sf.fastGenes.aa.fasta
        cat ./slowGenes/*.fa | grep -A 1 Sfur | grep -v "^\-\-" > sf.slowGenes.aa.fasta
        
        tblastn -db sf_GCA_014356515.1_genomic.fa -query sf.fastGenes.aa.fasta -out sf.fastGenes.aa.fasta.blast -evalue 1e-5 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 8
        tblastn -db sf_GCA_014356515.1_genomic.fa -query sf.slowGenes.aa.fasta -out sf.slowGenes.aa.fasta.blast -evalue 1e-5 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 8
        
        less -S sf.fastGenes.aa.fasta.blast | awk '($5 >= 80){print;}' > sf.fastGenes.aa.fasta.blast.filter
        less -S sf.slowGenes.aa.fasta.blast | awk '($5 >= 80){print;}' > sf.slowGenes.aa.fasta.blast.filter
        
        

## 2 - Distribution of ETLVEs in planthopper individuals

### 2.1 - Mapping

To understand the distribution of ETLVEs, we mapped Illumina reads to planthopper genomes using BWA MEM (version 0.7.17). Note that for your analysis, please replace $sampleID with the actual sample names. 

        # BWA MEM version 0.7.17
        ~/bin/bwa mem -M -t 40 -R '@RG\tID:'"$sampleID"'\tSM:'"$sampleID"'\tPL:Illumina' -o $sampleID.sam LsHau_GCA_014465815.1_genomic.fa ./links/$sampleID_1_clean.fq.gz ./links/$sampleID_2_clean.fq.gz
        
        # Convert SAM files to BAM files with SAMTools version 1.6
        ~/bin/samtools sort --threads 40 $sampleID.bam > $sampleID.sorted.bam
        ~/bin/samtools index -@ 40 $sampleID.sorted.bam
        
        # MarkDuplications using GATK version 4.2.6.1
        ~/bin/gatk-4.2.6.1/gatk MarkDuplicates -I $sampleID.sorted.bam -O $sampleID_markdup.bam -M $sampleID_markdup_metrics.txt
        
        # Index BAM files with SAMTools version 1.6
        ~/bin/samtools sort --threads 40 $sampleID_markdup.bam > $sampleID_markdup.sorted.bam
        ~/bin/samtools index -@ 40 $sampleID_markdup.sorted.bam
        
### 2.2 - Estimating overall sample depth with 100-bp sliding window using mosdepth

Mosdepth version 0.3.3 was used to estimate average sequencing depth for each sample. Note that please replace $sampleID with the actual sample names.

        for i in sf_GCA_014356515.1_genomic.fa ls_GCA_014465815.1_genomic.fa nl_GCA_014356525.1_genomic.fa
        do
         ## make a genome file
         samtools faidx $i
         awk -v OFS='\t' {'print $1,$2'} $i\.fai > genome.txt
        
         ## make 100 bp sliding window
         bedtools makewindows -g genome.txt -w 100 > $i\.bed.100
        
         # step 2 calculating depth using mosdepth
         # changed mapq from 20 to 0
         mosdepth --threads 20 --fast-mode --mapq 0 -n --by $i\.bed.100 $sampleID\.100bp $sampleID\.sorted.bam
        done

### 2.3 - Estimating per-base pair depth using mosdepth

To evaluate whether the EVEs are present in each sample, we used mosdepth to estimate per-base pair sequencing depth for each ETLVE. The ETLVE was identified as presence if the sequencing depth is higher than 5 for at least 50% of that ETLVE region. 
        
        # make 1 bp sliding window
        bedtools makewindows -g genome.txt -w 1 > SfHau_GCA_014356515.1_genomic.fa.bed.1
        
        # Prepare a bed file for each ETLVE, EVEs_sf.blast was from step1.1
        # For example, for Sogatella furcifera, the BLAST input file (generated from step1.1) looks like this: 
        # 2944	41453656	SfETLVE1_Chr7_length_2944_32054052_32056995	CM025296.1	100.000	2944	0	0	1	2944	32054052	32056995	0.0	5437
        # 1020	33928560	SfETLVE2_Chr8_length_1020_17350607_17351626	CM025297.1	100.000	1020	0	0	1	1020	17350607	17351626	0.0	1884
        # 822	33928560	SfETLVE3_Chr8_length_822_21854812_21855633	CM025297.1	100.000	822	0	0	1	822	21854812	21855633	0.0	1519
        # 822	34689497	SfETLVE3_Chr8_length_822_21854812_21855633	CM025300.1	99.859	709	1	0	1	709	16941783	16942491	0.0	1304
        
        # Create bed files
        perl batch_prep_mosdepth_bed.pl EVEs_sf.blast sf_GCA_014356515.1_genomic.fa.bed.1 sf_GCA_014356515.1_genomic.fa.bed.1
        
        mosdepth --threads 20 -b sf_GCA_014356515.1_genomic.fa.bed.1.EVE1 --fast-mode --mapq 0 $sampleID $sampleID\.sorted.bam
        mosdepth --threads 20 -b sf_GCA_014356515.1_genomic.fa.bed.1.EVE2 --fast-mode --mapq 0 $sampleID $sampleID\.sorted.bam
        mosdepth --threads 20 -b sf_GCA_014356515.1_genomic.fa.bed.1.EVE3 --fast-mode --mapq 0 $sampleID $sampleID\.sorted.bam
        mosdepth --threads 20 -b sf_GCA_014356515.1_genomic.fa.bed.1.EVE4 --fast-mode --mapq 0 $sampleID $sampleID\.sorted.bam
        
        # Calculate the percentage of ETLVE regions with depth >= 5x and coverage >= 50%
        for j in {1..4}
        do
         perl EVE_filter_depth.pl ./EVE$j\/$i\.regions.bed.gz ./EVE$j\/$i\.regions.bed.EVEfilter
        done
        
### 2.2 - SNP calling using GATK haplotypeCaller and create VCF files
        
        # HaplotypeCaller to make GVCF files
        ~/bin/gatk-4.2.6.1/gatk HaplotypeCaller -R NlHau_GCA_014356525.1_genomic.fa -I $sampleID_markdup.sorted.bam -O $sampleID_markdup.EVE.gvcf.gz --native-pair-hmm-threads 40 --intervals EVE_fast_slow.sort.list -ERC GVCF --disable-read-filter MappingQualityReadFilter
        
        # Combine GVCF files
        find . -type f -name "*_markdup.EVE.gvcf.gz" > input.gvcf.lis
        ~/bin/gatk-4.2.6.1/gatk CombineGVCFs -R LsHau_GCA_014465815.1_genomic.fa --variant input.gvcf.list -O cohort.EVE.g.vcf.gz
        
        # Genotype GVCF files
        ~/bin/gatk-4.2.6.1/gatk GenotypeGVCFs -R /share/home/liyiyuan/projects/collab_0_EVE/sample3_Laod/LsHau_GCA_014465815.1_genomic.fa -V cohort.EVE.g.vcf.gz -O output.EVE.g.vcf.gz
        
        # Filter VCF files
        ~/bin/gatk-4.2.6.1/gatk VariantFiltration -V output.EVE.g.vcf.gz -filter "AF < 0.1" --filter-name "AF0.1" -O output.EVE.snp_indel_filtered.vcf.gz
        
### 2.3 - Perl script to calculate heterozygosity
        # For each EVE
        # Example BLAST input: 
        # 2944	41453656	SfETLVE1_Chr7_length_2944_32054052_32056995	CM025296.1	100.000	2944	0	0	1	2944	32054052	32056995	0.0	5437
        perl LoP_calculation.all_gene_combined_hetero.v0.9.pl EVEs_sf.blast.modified.EVE1 all_chr_combined.vcf EVEs_sf.blast.modified.EVE1.LoP
        perl pick_LoP_on_list.pl EVEs_sf.blast.modified.EVE1.LoP EVE1.summary 0.5 sup4_EVEs_sf.blast.modified.EVE1.LoP.exist_samples
         done
