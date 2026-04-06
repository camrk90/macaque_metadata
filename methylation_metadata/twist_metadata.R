library(tidyverse)

#Generate functions-------------------------------------------------------------
edit_dnam<- function(x) {
  
  dd<- subset(x, select=twist_metadata$vantage_id)
  
  rownames(dd)<- gsub("Region_", "", rownames(dd))
  chrs<- str_split_i(rownames(dd), "\\.", 1)

  dd<- as.data.frame(dd) 
  dd$chr<- chrs
  
  dd_list<- split(dd, dd$chr)
  
  dd_list<- dd_list[-1]
  
  names(dd_list)<- 1:21
  
  dd_list <- lapply(dd_list, function(x) {
    select(x, -chr)
  })
  
  return(dd_list)
  
}

#Import metadata----------------------------------------------------------------
twist_metadata<- read_tsv("/scratch/sbaptis7/PrimateLifespan/rhesus/macaquelifespan_tms_metadata.tsv")
bad_sex<- read_csv("/scratch/sbaptis7/PrimateLifespan/VantageGenomic/inconsistent_sex_vantagesamples.csv")
bad_conv<- read_csv("/scratch/sbaptis7/PrimateLifespan/VantageGenomic/vantage_id_faulty_conversion.csv")

#Import methylation data--------------------------------------------------------
twist<- readRDS("/scratch/sbaptis7/PrimateLifespan/rhesus/rhesus_primlife_filtered_regions.rds")
meth<- twist[['methylation']]
cov<- twist[['coverage']]
regions<- twist[["regions"]]

chrs<- data.frame(chr=unique(regions$chr))
chrs$chrs<- c(1:22)
regions2<- left_join(regions, chrs, by='chr')
regions2<- regions2 %>% mutate(chr = chrs) %>% select(-chrs)
regions2<- regions2 %>% select(chr, start, end)
write.table(regions2, "/home/ckelsey4/Cayo_meth/twist_regions_short.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Remove samples with inconsistent sex and id assignments
bad_ids<- c(bad_sex$library, bad_conv$vantage_id)

twist_metadata<- twist_metadata[!twist_metadata$vantage_id %in% bad_ids,]

twist_metadata<- twist_metadata[twist_metadata$vantage_id %in% colnames(cov),]
twist_metadata$subject_id<- gsub("cayo:", "", twist_metadata$subject_id)
twist_metadata$tissue_id<- gsub("smacklab:", "", twist_metadata$tissue_id)

#Add n, within_age, and mean_age cols
twist_metadata<- twist_metadata %>%
  group_by(subject_id) %>%
  mutate(n=n()) %>%
  filter(n > 1) %>%
  mutate(age = as.numeric((capture_date - subject_dob)/365.25),
         mean_age = mean(age),
         within_age = age - mean_age) 

#Move important cols to front
twist_metadata<- twist_metadata %>%
  relocate(subject_id, subject_sex, age, within_age, mean_age, subject_dob, n)

#Subset cov/meth with metadata and make into chromosome list
cov_list<- edit_dnam(cov)
meth_list<- edit_dnam(meth)

#Write output files
write.table(twist_metadata, "/home/ckelsey4/Cayo_meth/twist_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(cov_list, "/home/ckelsey4/Cayo_meth/twist_cov_list")
saveRDS(meth_list, "/home/ckelsey4/Cayo_meth/twist_meth_list")

#Plot metadata 
twist_metadata<- twist_metadata %>%
  group_by(subject_id) %>%
  mutate(min_age = min(age))
twist_metadata$age<- round(twist_metadata$age, 0)

twist_metadata %>%
  ggplot(aes(x=age, y=reorder(subject_id, min_age), colour=as.factor(subject_sex))) +
  geom_path(linewidth = 0.5) +
  geom_point(colour="black", size = 0.25) +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_colour_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size=12) +
  theme(legend.position = "none", 
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(1, 1, 1, 1, "pt"))












