library(tidyverse)
library(GENESIS)

path<- "/scratch/ckelsey4/Cayo_meth/"

#Import metadata
blood_metadata<- read.table(paste(path, "blood_metadata_full.txt", sep=""))

rna_meta<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_metadata_wELA_18Jul25.rds")

rna_meta<- rna_meta %>%
  filter(Stimulation == "H2O") %>%
  select(c(animal_ID, sex, trapped_age, trapping_ID, Sample_ID, Seq_batch)) %>%
  arrange(Sample_ID)

#Generate vectors of animal biosample ids
dnam_animal_ids<- blood_metadata$monkey_id
dnam_bio_ids<- blood_metadata$lid_pid

rna_animal_ids<- rna_meta$animal_ID

#Import king output file
file.king <- c(paste(path, "king.kin0", sep=""))

generate_kinship<- function(animal_ids, bio_ids, bio_type){
  
  #Generate kin matrix
  kin.matrix<- kingToMatrix(file.king, estimator = "Kinship", sample.include=animal_ids)
  kinmat<- as.matrix(kin.matrix)
  
  #Arrange kinmat colnames by metadata to match r_matrix
  kinmat<- kinmat[unique(animal_ids),unique(animal_ids)]
  
  if (bio_type == "rna"){
    
    #generate empty z-matrix
    r_matrix<- data.frame(matrix(ncol = length(unique(animal_ids)), nrow=length(unique(animal_ids))))
    colnames(r_matrix)<- unique(animal_ids)
    rownames(r_matrix)<- unique(animal_ids)
    
    #Assign 1's to colnames that match ids in the id column (i.e. 1's for the same id)
    r_matrix[sapply(colnames(r_matrix), `==`, rownames(r_matrix))] <- 1
    
    r_matrix<- as.matrix(r_matrix)
    
    #Replace NAs with 0s
    r_matrix[is.na(r_matrix)]<- 0
    
    #Multiply matrices together to get full kinship matrix
    full_kin<- r_matrix %*% kinmat #%*% t(r_matrix)
    
  } else if (bio_type == "dnam") {
    
    #generate empty z-matrix
    r_matrix<- data.frame(matrix(ncol = length(unique(animal_ids)), nrow=length(bio_ids)))
    colnames(r_matrix)<- unique(animal_ids)
    rownames(r_matrix)<- bio_ids
    
    #Add monkey_id column to match colnames to
    all.equal(rownames(r_matrix), bio_ids)
    r_matrix$ids<- animal_ids
    
    #Assign 1's to colnames that match ids in the id column (i.e. 1's for the same id)
    r_matrix[sapply(colnames(r_matrix), `==`, r_matrix$ids)] <- 1
    
    #Remove id column and assign class matrix
    r_matrix<- r_matrix[,-length(r_matrix)]
    
    r_matrix<- as.matrix(r_matrix)
    
    #Replace NAs with 0s
    r_matrix[is.na(r_matrix)]<- 0
    
    #Multiply matrices together to get full kinship matrix
    full_kin<- r_matrix %*% kinmat %*% t(r_matrix)
    
  }
  
  return(full_kin)
  
}

dnam_kin<- generate_kinship(dnam_animal_ids, dnam_bio_ids)
rna_kin<- generate_kinship(rna_animal_ids, rna_bio_ids, "rna")


#Save output file as rds
saveRDS(dnam_kin, paste0(path, "full_kin_matrix.rds", sep=""))
saveRDS(rna_kin, paste(path, "rna_kin_matrix.rds", sep=""))
