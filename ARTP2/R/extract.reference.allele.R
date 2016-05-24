
extract.reference.allele <- function(stat){
  
  msg <- paste("Extracting allele information:", date())
  message(msg)
  
  snps <- NULL
  nstudy <- length(stat)
  for(i in 1:nstudy){
    snps <- unique(c(snps, stat[[i]][, 'SNP']))
  }
  
  nsnp <- length(snps)
  
  RefAllele <- rep(NA, nsnp)
  EffectAllele <- rep(NA, nsnp)
  names(RefAllele) <- snps
  names(EffectAllele) <- snps
  
  for(i in 1:nstudy){
    s <- stat[[i]][, 'SNP']
    RefAllele[s] <- stat[[i]][, 'RefAllele']
    EffectAllele[s] <- stat[[i]][, 'EffectAllele']
  }
  
  list(RefAllele = RefAllele, EffectAllele = EffectAllele)
  
}
