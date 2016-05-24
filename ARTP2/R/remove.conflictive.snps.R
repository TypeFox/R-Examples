
remove.conflictive.snps <- function(stat, ref.allele, conf.snps){
  
  if(!is.null(conf.snps)){
    nstudy <- length(stat)
    for(i in 1:nstudy){
      stat[[i]] <- stat[[i]][!(stat[[i]][, 'SNP'] %in% conf.snps), ]
    }
    
    RefAllele <- ref.allele$RefAllele
    ref.allele$RefAllele <- RefAllele[!(names(RefAllele) %in% conf.snps)]
    EffectAllele <- ref.allele$EffectAllele
    ref.allele$EffectAllele <- EffectAllele[!(names(EffectAllele) %in% conf.snps)]
  }
  
  list(stat = stat, ref.allele = ref.allele)
  
}
