
create.super.pathway.setup <- function(sum.stat, allele.info, reference, sub.pathway, options){
  
  nsnps <- nrow(allele.info)
  b <- (nsnps < 300000) # assume 1000 samples in reference data, then size of ref.geno is about 1.2 GB
  
  if(b){
    ref.geno <- load.reference.geno(reference, allele.info$SNP, options)
  }else{
    ref.geno <- NULL
  }
  
  setups <- list()
  ngrp <- length(sum.stat)
  # recover the summary statistics
  for(i in 1:ngrp){
    snps <- sum.stat[[i]]$snps.in.study
    ai <- allele.info[allele.info$SNP %in% snps, ]
    
    if(b){
      snp.id <- which(colnames(ref.geno) %in% snps)
      norm.stat <- recover.stat(sum.stat[[i]], sub.pathway[[i]], ref.geno[, snps, drop = FALSE], ai, options)
      ref.geno <- ref.geno[, -snp.id, drop = FALSE]
      gc()
    }else{
      ref.id <- unique(ai$Reference.ID)
      ref.geno <- load.reference.geno(reference[ref.id, ], sub.pathway[[i]]$SNP, options)
      norm.stat <- recover.stat(sum.stat[[i]], sub.pathway[[i]], ref.geno, ai, options)
      rm(ref.id, ref.geno)
      gc()
    }
    
    setups[[i]] <- list(options = options, allele.info = ai, pathway = sub.pathway[[i]], norm.stat = norm.stat)
    rm(snps, ai, norm.stat)
    gc()
  }
  
  setup <- merge.setups(setups)
  rm(setups)
  gc()
  
  setup
  
}

