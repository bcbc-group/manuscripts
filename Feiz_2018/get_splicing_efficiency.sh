#!/bin/bash

#Collect variables from script                
bam_name=$1;
splice=$2;
exon=$3;
intron=$4;
prefix=$5;

# strategy: run bedtools intersect to extract all the reads intersecting the exon/intron junction, in 3 steps: 
#1 find all the reads that intersect the splice sites, use -s to force strandeness
#2 eliminate purely exonic reads 
#3 eliminate purely intronic reads

bedtools intersect -s -a $bam_name -b $splice | bedtools intersect -a stdin -b $intron | bedtools intersect -a stdin -b $exon > $prefix.splice_junctions.bam
#bedtools intersect -s -split -a $bam_name -b $splice > $prefix.splice_junctions.bam 

#use samtools to only extract the reads that are spliced (exon/exon) from the splice_junctions.bam file and calculate the coverage using bedtools coverage
samtools view -h $prefix.splice_junctions.bam | awk '$6 ~ /N/ || $1 ~ /^@/' | samtools view -bS - | bedtools coverage -bed -a $intron -b stdin |cut -f 1,2,3,4,5,10 > $prefix.spliced_coverage.txt
#samtools view -h $prefix.splice_junctions.bam | awk '$6 ~ /N/ || $1 ~ /^@/' |samtools view -bS - | bedtools intersect -a $exon -b stdin | bedtools coverage -bed -a $splice -b stdin |cut -f 1,2,3,4,5,10 > $prefix.spliced_coverage.txt
#this should be changed to check the exact location of the splicing
#samtools view -h $prefix.splice_junctions.bam | awk '$6 ~ /N/ || $1 ~ /^@/' |samtools view -bS - | bedtools coverage -bed -a $splice -b stdin |cut -f 1,2,3,4,5,10 > $prefix.spliced_coverage.txt 

#use samtools to only extract the reads that are NOT spliced (exon/intron) from the splice_junctions.bam file and calculate the coverage using bedtools coverage
samtools view -h $prefix.splice_junctions.bam | awk '$6 !~ /N/ || $1 ~ /^@/' | samtools view -bS - | bedtools coverage -bed -a $intron -b stdin |cut -f 1,2,3,4,5,10 > $prefix.unspliced_coverage.txt
#samtools view -h $prefix.splice_junctions.bam | awk '$6 !~ /N/ || $1 ~ /^@/' | samtools view -bS - |bedtools intersect -a $intron -b stdin| bedtools coverage -bed -a $splice -b stdin |cut -f 1,2,3,4,5,10 > $prefix.unspliced_coverage.txt
#samtools view -h $prefix.splice_junctions.bam | awk '$6 !~ /N/ || $1 ~ /^@/' | samtools view -bS - | bedtools coverage -bed -a $splice -b stdin |cut -f 1,2,3,4,5,10 > $prefix.unspliced_coverage.txt
