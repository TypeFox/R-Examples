
update.ref.geno <- function(ref.geno, exc.snps){
  
  if(length(exc.snps) == 0){
    return(ref.geno)
  }
  
  id <- which(!(colnames(ref.geno) %in% exc.snps))
  if(length(id) == 0){
    msg <- "All SNPs excluded in update.ref.geno"
    stop(msg)
  }
  
  ref.geno <- ref.geno[, id, drop = FALSE]
  ref.geno
  
}

