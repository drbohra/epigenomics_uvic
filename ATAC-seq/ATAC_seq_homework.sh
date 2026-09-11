#!/bin/bash
set -e

mkdir -p data/bigBed.files data/bed.files analyses/peaks.analysis annotation

# Protein-coding promoter and gene-body coordinates
cp -f ../ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed annotation/
cp -f ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed annotation/

# Select ATAC-seq pseudoreplicated narrowPeak files from ENCODE metadata
grep -F "bigBed_narrowPeak" metadata.tsv | \
grep -F "pseudoreplicated_peaks" | \
grep -F "GRCh38" | \
awk 'BEGIN{FS=OFS="\t"}{print $1,$11,$5}' | \
sort -k2,2 -k1,1r | sort -k2,2 -u \
> analyses/bigBed.peaks.ids.txt

# Download peaks
cut -f1 analyses/bigBed.peaks.ids.txt | while read filename; do
    wget -O data/bigBed.files/"$filename".bigBed \
    "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

# Convert bigBed to BED
cut -f1 analyses/bigBed.peaks.ids.txt | while read filename; do
    bigBedToBed data/bigBed.files/"$filename".bigBed \
    data/bed.files/"$filename".bed
done

# Peaks overlapping promoters
cut -f1,2 analyses/bigBed.peaks.ids.txt | while read filename tissue; do
    bedtools intersect \
    -a data/bed.files/"$filename".bed \
    -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
    -u > analyses/peaks.analysis/"$tissue".ATAC.promoter.peaks.bed
done

# Peaks outside protein-coding gene bodies
cut -f1,2 analyses/bigBed.peaks.ids.txt | while read filename tissue; do
    bedtools intersect \
    -a data/bed.files/"$filename".bed \
    -b annotation/gencode.v24.protein.coding.gene.body.bed \
    -v > analyses/peaks.analysis/"$tissue".ATAC.outside.genes.bed
done

echo "Promoter-overlapping peaks:"
wc -l analyses/peaks.analysis/*.ATAC.promoter.peaks.bed

echo "Peaks outside genes:"
wc -l analyses/peaks.analysis/*.ATAC.outside.genes.bed
