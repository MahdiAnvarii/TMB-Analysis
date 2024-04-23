#!/bin/bash

# Step 1 : QC - Run Fastqc
fastqc /mnt/d/NGS/Samples/P2/Reads/SRR28000175_1.fastq -o /mnt/d/NGS/Samples/P2/Reads/
fastqc /mnt/d/NGS/Samples/P2/Reads/SRR28000175_2.fastq -o /mnt/d/NGS/Samples/P2/Reads/

# Step 2 : Trimming - Run Trimmomatic
trimmomatic PE -phred33 /mnt/d/NGS/Samples/P2/Reads/SRR28000175_1.fastq /mnt/d/NGS/Samples/P2/Reads/SRR28000175_2.fastq /mnt/d/NGS/Samples/P2/Reads/SRR28000175_forward_paired.fastq.gz /mnt/d/NGS/Samples/P2/Reads/SRR28000175_forward_unpaired.fastq.gz /mnt/d/NGS/Samples/P2/Reads/SRR28000175_reverse_paired.fastq.gz /mnt/d/NGS/Samples/P2/Reads/SRR28000175_reverse_unpaired.fastq.gz ILLUMINACLIP:/home/mahdi/mydir/Packages/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

fastqc /mnt/d/NGS/Samples/P2/Reads/SRR28000175_forward_paired.fastq.gz -o /mnt/d/NGS/Samples/P2/Reads/
fastqc /mnt/d/NGS/Samples/P2/Reads/SRR28000175_reverse_paired.fastq.gz -o /mnt/d/NGS/Samples/P2/Reads/

# Step 3 : Map to reference - Use BWA-MEM
bwa mem -t 6 -R "@RG\tID:SRR2800017\tPL:ILLUMINA\tSM:SRR2800017" /mnt/d/NGS/References/hg38.fa /mnt/d/NGS/Samples/P2/Reads/SRR28000175_forward_paired.fastq.gz /mnt/d/NGS/Samples/P2/Reads/SRR28000175_reverse_paired.fastq.gz > /mnt/d/NGS/Samples/P2/Aligned/SRR28000175.paired.sam
samtools view /mnt/d/NGS/Samples/P2/Aligned/SRR28000175.paired.sam | less
samtools flagstat /mnt/d/NGS/Samples/P2/Aligned/SRR28000175.paired.sam

# Step 4 : Mark Duplicates and Sort - Run GATK4
gatk MarkDuplicatesSpark -I /mnt/d/NGS/Samples/P2/Aligned/SRR28000175.paired.sam -O /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup.bam
samtools flagstat /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup.bam

# Step 5 : Base Quality Score Recalibration
gatk BaseRecalibrator -I /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup.bam -R /mnt/d/NGS/References/hg38.fa --known-sites /mnt/d/NGS/References/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf -O /mnt/d/NGS/Samples/P2/Aligned/recal_data.table
gatk ApplyBQSR -I /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup.bam -R /mnt/d/NGS/References/hg38.fa --bqsr-recal-file /mnt/d/NGS/Samples/P2/Aligned/recal_data.table -O /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup_bqsr.bam

# Step 6 : Collect Alignment and Insert Size Metrics
gatk CollectAlignmentSummaryMetrics R=/mnt/d/NGS/References/hg38.fa I=/mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup_bqsr.bam O=/mnt/d/NGS/Samples/P2/Aligned/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=/mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup_bqsr.bam OUTPUT=/mnt/d/NGS/Samples/P2/Aligned/insert_size_metrics.txt HISTOGRAM_FILE=/mnt/d/NGS/Samples/P2/Aligned/insert_size_histogram.pdf
multiqc .

# Step 7 : Call Variants - Run GATK4 Haplotype Caller
gatk HaplotypeCaller -R /mnt/d/NGS/References/hg38.fa -I /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup_bqsr.bam -O /mnt/d/NGS/Samples/P2/Results/raw_variants.vcf

# Step 8 : Extract SNPs and INDELs
gatk SelectVariants -R /mnt/d/NGS/References/hg38.fa -V /mnt/d/NGS/Samples/P2/Results/raw_variants.vcf --select-type SNP -O /mnt/d/NGS/Samples/P2/Results/raw_snps.vcf
gatk SelectVariants -R /mnt/d/NGS/References/hg38.fa -V /mnt/d/NGS/Samples/P2/Results/raw_variants.vcf --select-type INDEL -O /mnt/d/NGS/Samples/P2/Results/raw_indels.vcf

# Step 9 : Filter Variants - GATK4
gatk VariantFiltration -R /mnt/d/NGS/References/hg38.fa -V /mnt/d/NGS/Samples/P2/Results/raw_snps.vcf -O /mnt/d/NGS/Samples/P2/Results/filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"
gatk VariantFiltration -R /mnt/d/NGS/References/hg38.fa -V /mnt/d/NGS/Samples/P2/Results/raw_indels.vcf -O /mnt/d/NGS/Samples/P2/Results/filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0" -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

# Step 10 : Select Variants that pass filters
gatk SelectVariants --exclude-filtered -V /mnt/d/NGS/Samples/P2/Results/filtered_snps.vcf -O /mnt/d/NGS/Samples/P2/Results/analysis-ready-snps.vcf
gatk SelectVariants --exclude-filtered -V /mnt/d/NGS/Samples/P2/Results/filtered_indels.vcf -O /mnt/d/NGS/Samples/P2/Results/analysis-ready-indels.vcf

cat /mnt/d/NGS/Samples/P2/Results/analysis-ready-snps.vcf | grep -v -E "DP_filter|GQ_filter" > /mnt/d/NGS/Samples/P2/Results/analysis-ready-snps-GT.vcf
cat /mnt/d/NGS/Samples/P2/Results/analysis-ready-indels.vcf | grep -v -E "DP_filter|GQ_filter" > /mnt/d/NGS/Samples/P2/Results/analysis-ready-indels-GT.vcf

# Step 11 : Annotate Variants - GATK Funcotator
gatk Funcotator --variant /mnt/d/NGS/Samples/P2/Results/analysis-ready-snps-GT.vcf --reference /mnt/d/NGS/References/hg38.fa --ref-version hg38 --data-sources-path /mnt/d/NGS/References/funcotator_dataSources.v1.7.20200521g/ --output /mnt/d/NGS/Samples/P2/Results/analysis-ready-snps-GT-funcotated.vcf --output-file-format VCF
gatk Funcotator --variant /mnt/d/NGS/Samples/P2/Results/analysis-ready-indels-GT.vcf --reference /mnt/d/NGS/References/hg38.fa --ref-version hg38 --data-sources-path /mnt/d/NGS/References/funcotator_dataSources.v1.7.20200521g/ --output /mnt/d/NGS/Samples/P2/Results/analysis-ready-indels-GT-funcotated.vcf --output-file-format VCF

# Step 12 : TMB Calculation
bedtools bamtobed -i /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup_bqsr.bam > /mnt/d/NGS/Samples/P2/Data/SRR28000175_sorted_dedup_bqsr.bed
bedtools merge -i /mnt/d/NGS/Samples/P2/Data/SRR28000175_sorted_dedup_bqsr.bed > /mnt/d/NGS/Samples/P2/Data/SRR28000175_sorted_dedup_bqsr_merged.bed
indels_count=$(grep -vc '^#' /mnt/d/NGS/Samples/P2/Results/analysis-ready-indels-GT-funcotated.vcf)
snps_count=$(grep -vc '^#' /mnt/d/NGS/Samples/P2/Results/analysis-ready-snps-GT-funcotated.vcf)
total_mutations=$((indels_count + snps_count))
# echo $total_mutations
python TMB.py $total_mutations

# Step 13 : MSI CAlculation
msisensor scan -d /mnt/d/NGS//References/hg38.fa -o /mnt/d/NGS/References/MSIscan.bed
msisensor msi -d /mnt/d/NGS/References/MSIscan.bed -t /mnt/d/NGS/Samples/P2/Aligned/SRR28000175_sorted_dedup_bqsr.bam -o /mnt/d/NGS/Samples/P2/Data/MSIoutput.bed
