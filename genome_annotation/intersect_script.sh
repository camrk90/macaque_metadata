#!/bin/sh

#This script intersects CpG sites with pre-defined genomic regions and then with an annotation file.
#Annotation files should be in .bed format and located in the same directory as this script.

module load bedtools2-2.30.0-gcc-11.2.0

#Input files
CPGs=$1
REGIONS=$3

#Output files
CPGSxREGIONS=regions_to_cpgs.txt

#Intersect autosome/x-chrom cpgs and regions
bedtools intersect -a ${CPGs} -b ${REGIONS} -wa -wb > ${CPGSxREGIONS}

#intersect with the annotation file
for file in *.bed; do
    echo "Intersecting ${file} with annotation file..."

    bedtools intersect -a ${file} -b ${CPGSxREGIONS} -wo > ${file%.bed}_annotation_intersect.txt
done

