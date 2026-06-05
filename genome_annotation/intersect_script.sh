#!/bin/sh

#This script intersects CpG sites with pre-defined genomic regions and then with an annotation file.
#Annotation files should be in .bed format and located in the same directory as this script.

module load bedtools2-2.30.0-gcc-11.2.0

#Set file paths
out_loc="$HOME/Cayo_meth/intersect_files/"
inputs_loc="$HOME/Projects/macaque_metadata/output_files/"

#Input files
CPGs=${inputs_loc}$1
REGIONS=${inputs_loc}$2

#Output files
CPGSxREGIONS=${out_loc}regions_to_cpgs.txt

#Intersect autosome/x-chrom cpgs and regions
bedtools intersect -a ${CPGs} -b ${REGIONS} -wa -wb > ${CPGSxREGIONS}

#intersect with the annotation file
for file in ${inputs_loc}mmul_*.bed; do
    base=$(basename "${file}" .bed)

    echo "Intersecting ${file} with annotation file..."

    bedtools intersect \
     -a ${file} -b ${CPGSxREGIONS} \
     -wo > ${out_loc}${base}_intersect.txt
done

