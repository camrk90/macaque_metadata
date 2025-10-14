#Gets M/Cov from regions, promoter, gene, or any other coordinates specified
dnam_to_regions<- function(bsseq_list, region_list, coverage_type, chr_names){
  
  return_list<- parallel::mclapply(names(bsseq_list),function(x){
    dd=getCoverage(cayo_filtered_list[[x]], regions = region_list[region_list@seqnames==x,], type = coverage_type, 
                   what = "perRegionTotal")
    rownames(dd)=region_list[region_list@seqnames==x,]$gene_id
    return(as.data.frame(dd))
  },mc.cores=20)
  
  names(return_list)<- chr_names
  
  return(return_list)
}

#Filters for constitutively/non-constitutively methylated regions/proms/genes
filter_tails<- function(M, Cov, min_perc, max_perc){
  
  #Generate percent methylation matrix
  ratio<- mapply('/', M, Cov, SIMPLIFY = F)
  
  #Filter out methylated genes that are constituvively hyper/hypomethylated
  ratio<- sapply(names(ratio), function(i){
    
    df<- ratio[[i]]
    #Keep rows where the number of na's does not equal the number of columns
    df<- df[rowSums(is.na(df)) != ncol(df),]
    #Filter for %meth below 0.10 and above 0.90
    df<- df[!rowMeans(df, na.rm=T) > max_perc,] 
    df<- df[!rowMeans(df, na.rm=T) < min_perc,]
    
  },simplify = F, USE.NAMES = T)
  
  #Use filtered ratio df to filter coverage matrix
  cov_out<- sapply(names(Cov), function(i){
    
    cov_chr<- Cov[[i]]
    sub_chr<- ratio[[i]]
    cov_chr<- cov_chr[rownames(cov_chr) %in% rownames(sub_chr),]
    
  },simplify = F, USE.NAMES = T)
  
  return(cov_out)
}


#Filters samples with less than 70% of regions covered and regions with less than 25% coverage
filter_regions<- function(region_list){
  
  #Unlist chrs to full matrix of regions across the genome
  df<- do.call(rbind, region_list[c(1:20, "X")])
  
  #Remove NAs
  df<- df %>% drop_na()
  
  #Filter samples with less than 70% of regions covered and regions with less than 25% coverage
  df<- df[,colSums(df >= 1) >= 0.70*nrow(df)]
  df<- df[rowSums(df >= 1) >= 0.25*ncol(df),]
  
  #Generate col with chromosome number
  df$chr<- rownames(df)
  df$chr<- str_split_i(df$chr,"\\.",1)
  
  #Re-list regions
  r_list<- split(df[, -ncol(df)], df$chr)
  return(r_list)
}