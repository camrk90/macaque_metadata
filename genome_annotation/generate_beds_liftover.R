library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)

######################################
### Generate .bed file for regions ###
######################################
#Import metadata and regions rds------------------------------------------------
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")

cov<- do.call(rbind, regions_cov)
regions<- rownames(cov)
regions<- str_split_i(regions, "\\.", 4)
regions<- as.data.frame(regions) %>%
  separate_wider_delim(regions, "_", names = c("chrom", "chromStart", "chromEnd"))

#Add 'chr' string to chr col
regions$chrom<- gsub(" ", "", paste("chr", regions$chrom, ""))

#Generate .bed file
regions_bed<- GenomicRanges::makeGRangesFromDataFrame(regions)

######################################
###   Generate .bed file for TEs   ###
######################################
#Import repeatmasker repeat coordinates-----------------------------------------
repeats<- read.table("/scratch/ckelsey4/Cayo_meth/rmsk.txt", header = F)
repeats<- repeats[, c(6:8, 11:13)]
colnames(repeats)<- c("chrom", "chromStart", "chromEnd", "repName", "repClass", "repFamily")

#Removes unused chrs such as "chr10_NW_021160242v1_random", "chrUn", "chrY" & "chrM"
chrs<- unique(gsub("_.*", "", repeats$chrom))
chrs<- chrs[!chrs %in% c("chrUn", "chrY", "chrM")]
repeats<- repeats %>% filter(chrom %in% chrs)

#Add 1 to coordinates to account for the fact that exporting to bed will subtract 1 to make it 0-based
#the repeats file is already 0-based
repeats<- repeats %>%
  mutate(chromStart = chromStart + 1, chromEnd = chromEnd + 1)

#Generate range col to more easily match joins in analysis script
repeats<- repeats %>%
  mutate(range = paste(as.character(chromStart), "-", as.character(chromEnd)))

#Sort chrs numerically
repeats<- repeats[str_order(repeats$chrom, numeric = TRUE), ]

##########################################
### Generate .bed file for CHMM states ###
##########################################
#Import bed file----------------------------------------------------------------
bed_file<- ("/scratch/ckelsey4/Cayo_meth/E062_15_coreMarks_hg38lift_dense.bed")
chmm_bed<- rtracklayer::import(bed_file, format = "bed")

#Import ucsc chain file
hg38_mmul10<- rtracklayer::import.chain("/scratch/ckelsey4/Cayo_meth/hg38ToRheMac10.over.chain")

#Liftover hg38 coords to mmul10 assembly and sort
seqlevelsStyle(chmm_bed) = "UCSC"
chmm_mmul<- rtracklayer::liftOver(chmm_bed, hg38_mmul10)
chmm_mmul<- unlist(chmm_mmul)
genome(chmm_mmul) = "mmul_10"
chmm_mmul<- sortSeqlevels(chmm_mmul)
selected_chrs<- paste("chr", c(1:20, "X"), sep = "")
chmm_mmul2<- keepSeqlevels(chmm_mmul, selected_chrs, pruning.mode = "coarse")
chmm_mmul2<- sort(chmm_mmul2)

#Export files as .bed for intersect---------------------------------------------
rtracklayer::export.bed(regions_bed, con = "./genome_annotation/regions.bed")
rtracklayer::export.bed(repeats, con = "./genome_annotation/mmul_repeats.bed")
rtracklayer::export.bed(chmm_mmul2, con = "./genome_annotation/mmul_chmm.bed")
