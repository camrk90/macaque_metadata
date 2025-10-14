library(tidyverse)
library(GENESIS)

path<- "/scratch/ckelsey4/Cayo_meth/"

#Import metadata
blood_metadata<- read.table(paste(path, "blood_metadata_full.txt", sep=""))

#Generate vector of animal ids
monkey_vector<- blood_metadata$monkey_id

#Import king output file
file.king <- c("king.kin0")

#Generate kin matrix
kin.matrix<- kingToMatrix(file.king, estimator = "Kinship", sample.include=monkey_vector)
kinmat<- as.matrix(kin.matrix)

#Arrange kinmat colnames by metadata to match r_matrix
kinmat<- kinmat[unique(blood_metadata$monkey_id),unique(blood_metadata$monkey_id)]

#Generate id vectors for Z-matrix of samples(lids) x individuals(monkey_id)
monkey_ids<- blood_metadata %>%
  dplyr::select(monkey_id) %>%
  unique()
lids<- blood_metadata %>%
  ungroup() %>%
  dplyr::select(lid_pid)

#generate empty z-matrix
r_matrix<- data.frame(matrix(ncol = nrow(monkey_ids), nrow=nrow(lids)))
colnames(r_matrix)<- monkey_ids$monkey_id
rownames(r_matrix)<- lids$lid_pid

#Add monkey_id column to match colnames to
all.equal(rownames(r_matrix), blood_metadata$lid_pid)
r_matrix$ids<- blood_metadata$monkey_id

#Assign 1's to colnames that match ids in the id column (i.e. 1's for the same id)
r_matrix[sapply(colnames(r_matrix), `==`, r_matrix$ids)] <- 1

#Remove id column and assign class matrix
r_matrix<- r_matrix %>%
  dplyr::select(-length(r_matrix))
r_matrix<- as.matrix(r_matrix)

#Replace NAs with 0s
r_matrix[is.na(r_matrix)]<- 0

#Multiply matrices together to get full kinship matrix
full_kin<- r_matrix %*% kinmat %*% t(r_matrix)

#Save output file as rds
saveRDS(full_kin, paste0(path, "full_kin_matrix.rds", sep=""))
