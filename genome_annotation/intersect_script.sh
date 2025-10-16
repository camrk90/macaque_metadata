#!/bin/sh

module load bedtools2-2.30.0-gcc-11.2.0

#All input files should be in bed file format

#Input files
CPGs=$1
Y_CPGs=$2
REGIONS=$3
REGIONS_Y=$4
REPEATS=$5
CHMM=$6

#Output files
CPGSxREGIONS=regions_to_cpgs.txt
Y_CPGSxREGIONS=y_regions_to_cpgs.txt
ALL_CPGS=all_regions_to_cpgs.txt

#Intersect autosome/x-chrom cpgs and regions
bedtools intersect -a ${CPGs} -b ${REGIONS} -wa -wb > ${CPGSxREGIONS}

#Intersect y-chrom cpgs and regions
bedtools intersect -a ${Y_CPGs} -b ${REGIONS_Y} -wa -wb > ${Y_CPGSxREGIONS}

#Concatenate autosomes/x-chrom/y-chrom cpgs
cat ${CPGSxREGIONS} ${Y_CPGSxREGIONS} > ${ALL_CPGS}

#intersect with the repeats file
bedtools intersect -a ${REPEATS} -b ${ALL_CPGS} -wo > region_repeat_intersect.txt

#intersect with the lifted over chmm file
bedtools intersect -a ${CHMM} -b ${ALL_CPGS} -wo > region_chmm_intersect.txt
