
update.ref.snps <- function(ref.snps, exc.snps){
  
  if(length(exc.snps) == 0){
    return(ref.snps)
  }
  
  ref.snps <- setdiff(ref.snps, exc.snps)
  if(length(ref.snps) == 0){
    msg <- "All SNPs excluded in update.ref.snps"
    stop(msg)
  }
  ref.snps
  
}
