library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)

#Set chrs to keep
selected_chrs<- paste("chr", c(1:20, "X"), sep = "")

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
regions<- regions %>%
  mutate(chromStart = as.numeric(chromStart) + 1)

#Generate .bed file
regions_bed<- GenomicRanges::makeGRangesFromDataFrame(regions)

#################################################
###   Generate .bed file for Islands/Shores   ###
#################################################
#Load in CpG Islands file
cpg_islands<- read.csv("/home/ckelsey4/Cayo_meth/intersect_files/cpgIslandExt.txt", sep = "\t", header = F)
cpg_islands<- cpg_islands[, 2:4]
colnames(cpg_islands)<- c("chrom", "chromStart", "chromEnd")

#Remove bad chrom names
island_chrs<- unique(cpg_islands$chrom)
island_chrs<- island_chrs[!grepl("_NW_", island_chrs)]
cpg_islands<- cpg_islands %>%
  filter(chrom %in% island_chrs) %>%
  mutate(chromStart = chromStart + 1)

#Generate CpG shores 2kb upstream of island start site 
cpg_shores<- cpg_islands %>%
  dplyr::select(chrom, chromStart) %>% 
  dplyr::rename(chromEnd = chromStart) %>%
  mutate(chromStart = chromEnd - 2000) %>%
  relocate(chromStart, .before = chromEnd)

#Generate Cpg shores 2kb downstream of island end site
cpg_shores_1<- cpg_islands %>%
  dplyr::select(chrom, chromEnd) %>% 
  dplyr::rename(chromStart = chromEnd) %>%
  mutate(chromEnd = chromStart + 2000)

#Bind shores and make any negative values 1
cpg_shores<- rbind(cpg_shores, cpg_shores_1)
cpg_shores<- cpg_shores %>%
  mutate(chromStart = ifelse(chromStart < 0, 1, chromStart))
cpg_shores<- cpg_shores %>%
  arrange(chrom, chromStart)

##############################################
###   Generate .bed file for Genes/Proms   ###
##############################################
paths<- c("/scratch/ckelsey4/Cayo_meth/macaque_genes.rds", 
          "/scratch/ckelsey4/Cayo_meth/macaque_promoters.rds")

genes_proms<- lapply(setNames(paths, c("genes", "proms")), function(x){
  
  #Load in genes/proms files
  gr_obj<- readRDS(x)
  
  #Edit seqnames to include "chr"
  seqlevelsStyle(gr_obj)<- "UCSC"
  
  gr_obj<- keepSeqlevels(gr_obj, selected_chrs, pruning.mode = "coarse")
  
  gr_obj<- sort(gr_obj)
  
})

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

#Generate range col to more easily match joins in analysis script
repeats<- repeats %>%
  mutate(range = paste(as.character(chromStart), "-", as.character(chromEnd)))

#Add 1 to coordinates to account for the fact that exporting to bed will subtract 1 to make it 0-based
#the repeats file is already 0-based
repeats<- repeats %>%
  mutate(chromStart = chromStart + 1)

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
rtracklayer::export.bed(regions_bed, con = "./output_files/regions.bed")
rtracklayer::export.bed(repeats, con = "./output_files/mmul_repeats.bed")
rtracklayer::export.bed(chmm_mmul2, con = "./output_files/mmul_chmm.bed")
rtracklayer::export.bed(genes_proms[["genes"]], con = "./output_files/mmul_genes.bed")
rtracklayer::export.bed(genes_proms[["proms"]], con = "./output_files/mmul_proms.bed")
rtracklayer::export.bed(cpg_islands, con = "./output_files/mmul_islands.bed")
rtracklayer::export.bed(cpg_shores, con = "./output_files/mmul_shores.bed")

