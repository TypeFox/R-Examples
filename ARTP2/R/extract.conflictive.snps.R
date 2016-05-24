
extract.conflictive.snps <- function(stat, ref.allele){
  
  msg <- paste("Extracting SNPs with conflictive alleles:", date())
  message(msg)
  
  RefAllele <- ref.allele$RefAllele
  EffectAllele <- ref.allele$EffectAllele
  
  snps <- names(RefAllele)
  
  A1 <- paste(RefAllele, EffectAllele, sep = '')
  A2 <- paste(EffectAllele, RefAllele, sep = '')
  names(A1) <- snps
  names(A2) <- snps
  
  nstudy <- length(stat)
  conf.snps <- NULL
  for(i in 1:nstudy){
    a <- paste(stat[[i]][, 'RefAllele'], stat[[i]][, 'EffectAllele'], sep = '')
    s <- stat[[i]][, 'SNP']
    names(a) <- s
    b <- which(a != A1[s] & a != A2[s])
    if(length(b) == 0){
      next
    }
    conf.snps <- c(conf.snps, names(b))
  }
  
  conf.snps <- unique(conf.snps)
  conf.snps
  
  
}
