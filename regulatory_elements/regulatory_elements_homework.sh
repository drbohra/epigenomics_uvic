#!/bin/bash
set -e

# Task 5: distal regulatory element analysis

mkdir -p data/bigBed.files data/bed.files analyses

# Candidate distal regulatory elements:
# distal ATAC peaks overlapping both H3K27ac and H3K4me1

bedtools intersect \
-a ../ATAC-seq/analyses/peaks.analysis/sigmoid_colon.ATAC.outside.genes.bed \
-b data/bed.files/ENCFF872UHN.bed -u | \
bedtools intersect -a - -b data/bed.files/ENCFF724ZOF.bed -u \
> analyses/sigmoid_colon.regulatory.elements.bed

bedtools intersect \
-a ../ATAC-seq/analyses/peaks.analysis/stomach.ATAC.outside.genes.bed \
-b data/bed.files/ENCFF977LBD.bed -u | \
bedtools intersect -a - -b data/bed.files/ENCFF844XRN.bed -u \
> analyses/stomach.regulatory.elements.bed

# Chromosome 1 regulatory-element starts
cat analyses/sigmoid_colon.regulatory.elements.bed \
analyses/stomach.regulatory.elements.bed | \
awk 'BEGIN{FS=OFS="\t"} $1=="chr1"{print $4,$2}' \
> regulatory.elements.starts.tsv

# Chromosome 1 protein-coding gene transcription starts
awk 'BEGIN{FS=OFS="\t"} $1=="chr1"{
if($6=="+"){start=$2}else{start=$3};
print $4,start
}' ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed \
> gene.starts.tsv

# Closest gene for each regulatory element
while read element start; do
    result=$(python ../bin/get.distance.py \
        --input gene.starts.tsv \
        --start "$start")
    printf "%s\t%s\n" "$element" "$result"
done < regulatory.elements.starts.tsv \
> regulatoryElements.genes.distances.tsv

# Mean and median distance
Rscript -e 'x<-read.table("regulatoryElements.genes.distances.tsv",sep="\t",header=FALSE); cat("Mean distance:",mean(x[,4]),"\nMedian distance:",median(x[,4]),"\n")'
