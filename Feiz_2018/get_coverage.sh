#!/bin/bash

#Collect variables from script
name=$1
bam_name=$2
exon=$3
intron=$4
prefix=$5

#get counts
#counts="$(samtools view -F 256 $bam_name | wc -l)"
counts="$(samtools view -F 4 $bam_name | wc -l)"

echo "$counts" > $prefix.counts
echo "$bam_name" > $prefix.bam_name

#get strand specific, single base coverage using bedtools on the bam file from the tophat alignment. use -s for strandeness. First for plus strand, then minus strand. Normalized with the number of aligned reads. SS added split
bedtools genomecov -split -strand + -d -ibam $bam_name -g genome_ChrC.bed | grep "^$name" | awk -v awk_count=$counts '{print $1 "\t" $2 "\t" ($3*1000000)/(awk_count)}' > $prefix.coverage_plus.csv
bedtools genomecov -split -strand - -d -ibam $bam_name -g genome_ChrC.bed | grep "^$name" | awk -v awk_count=$counts '{print $1 "\t" $2 "\t" ($3*1000000)/(awk_count)}'> $prefix.coverage_minus.csv

#get strand specific coverage on a 100nt sliding window (every 50nt), normalized by the number of aligned reads.  SS added split
bedtools intersect -split -s -bed -wb -a $bam_name -b 100nt_plus.bed | cut -f 13,14,15 | sort -k2,2n | uniq -c |  awk -v awk_count=$counts '{print $2 "\t" $3 "\t" $4 "\t" ($1*1000000)/(awk_count)}' > $prefix.window_plus.csv
bedtools intersect -split -s -bed -wb -a $bam_name -b 100nt_minus.bed | cut -f 13,14,15 | sort -k2,2n | uniq -c | awk -v awk_count=$counts '{print $2 "\t" $3 "\t" $4 "\t" ($1*1000000)/(awk_count)}' > $prefix.window_minus.csv

#calculate RPKM for all the exons. SS added split
#bedtools intersect -split -s -bed -wb -a $bam_name -b $exon | cut -f 13,16,17 | sort -k2,2n | uniq -c |sed -r 's/^( *[^ ]+) +/\1\t/' > $prefix.exp_exon.txt
bedtools multicov -split -s -bams $bam_name -bed $exon |cut -f1,4,5,9,10 > $prefix.exp_exon.txt

#calculate RPKM for all the introns. SS added split
#bedtools intersect -split -s -bed -wb -a $bam_name -b $intron | cut -f 13,16,17 | sort -k2,2n | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > $prefix.exp_intron.txt
bedtools multicov -split -s -bams $bam_name -bed $intron |cut -f1,4,5,9,10 > $prefix.exp_intron.txt
