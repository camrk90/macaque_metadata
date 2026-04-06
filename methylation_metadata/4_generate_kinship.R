library(tidyverse)
library(GENESIS)

path<- "/home/ckelsey4/"

#Import metadata
blood_metadata<- read.table(paste(path, "Cayo_meth/twist_metadata.txt", sep=""), sep = "\t", header = T)
blood_metadata<- blood_metadata %>%
  arrange(vantage_id)

rna_meta<- read.table(paste(path, "rna_data/base_meta.txt", sep = ""))
rna_meta<- rna_meta %>%
  arrange(Sample_ID)

#Generate vectors of animal biosample ids
dnam_animal_ids<- blood_metadata$subject_id
dnam_bio_ids<- blood_metadata$vantage_id

rna_animal_ids<- rna_meta$animal_ID
rna_bio_ids<- rna_meta$Sample_ID

#Import king output file
file.king <- "king.kin0"

generate_kinship<- function(animal_ids, bio_ids, bio_type){
  
  #Generate kin matrix
  kin.matrix<- kingToMatrix("/scratch/ckelsey4/Cayo_meth/king.kin0", estimator = "Kinship", sample.include=animal_ids)
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

dnam_kin<- generate_kinship(dnam_animal_ids, dnam_bio_ids, "dnam")
rna_kin<- generate_kinship(rna_animal_ids, rna_bio_ids, "rna")

#Save output file as rds
saveRDS(dnam_kin, paste0(path, "Cayo_meth/dnam_kin_matrix.rds", sep=""))
saveRDS(rna_kin, paste(path, "rna_data/rna_kin_matrix.rds", sep=""))
