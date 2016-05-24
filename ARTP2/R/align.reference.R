
align.reference <- function(ref.geno, allele.info, options){
  
  msg <- paste("Realigning allele information of reference:", date())
  if(options$print) message(msg)
  
  rs <- allele.info$SNP
  ref.geno <- ref.geno[, rs, drop = FALSE]
  eaf <- apply(ref.geno, 2, function(u){mean(u, na.rm = TRUE)/2})
  id <- which(eaf > .5)
  
  if(length(id) > 0){
    ea <- allele.info$EffectAllele
    ra <- allele.info$RefAllele
    ea[id] <- allele.info$RefAllele[id]
    ra[id] <- allele.info$EffectAllele[id]
    allele.info$EffectAllele <- ea
    allele.info$RefAllele <- ra
  }
  
  foo <- function(g){
    f <- mean(g, na.rm = TRUE)/2
    if(f > .5){
      return(2 - g)
    }else{
      return(g)
    }
  }
  ref.geno <- apply(ref.geno, 2, foo)
  colnames(ref.geno) <- rs
  
  list(ref.geno = ref.geno, allele.info = allele.info)
  
}
