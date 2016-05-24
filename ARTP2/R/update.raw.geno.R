
update.raw.geno <- function(raw.geno, exc.snps){
  
  if(length(exc.snps) == 0){
    return(raw.geno)
  }
  
  id <- which(!(colnames(raw.geno) %in% exc.snps))
  if(length(id) == 0){
    msg <- "All SNPs excluded in update.raw.geno"
    stop(msg)
  }
  
  raw.geno <- raw.geno[, id, drop = FALSE]
  gc()
  raw.geno
  
}

