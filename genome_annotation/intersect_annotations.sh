#!/bin/sh

#module load bedtools2-2.30.0-gcc-11.2.0

ANNO=$1
ANNO_BED="${ANNO%%.*}.bed"
BED_FILE=$2
FILE_NAME="${BED_FILE%%.*}_intersect.txt"

#This chunk converts a .txt file to .bed file if need be and runs the intersect
if [[ ${ANNO} == *".txt" ]]; then
    cat ${ANNO} > ${ANNO_BED}
    #bedtools intersect -a ${BED_FILE} -b ${ANNO} -wo > ${FILE_NAME}
fi

#intersect with the lifted over chmm file
#bedtools intersect -a ${BED_FILE} -b ${ANNO} -wo > ${FILE_NAME}

