
update.sf <- function(sf, exc.snps){
  
  if(length(exc.snps) == 0){
    return(sf)
  }
  
  id <- which(!(sf$SNP %in% exc.snps))
  if(length(id) == 0){
    msg <- "All SNPs excluded in update.sf"
    stop(msg)
  }
  
  sf <- sf[!(sf$SNP %in% exc.snps), ]
  sf
  
}
