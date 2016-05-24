
update.allele.info <- function(allele.info, exc.snps){
  
  if(length(exc.snps) == 0){
    return(allele.info)
  }
  
  id <- which(!(allele.info$SNP %in% exc.snps))
  if(length(id) == 0){
    msg <- "All SNPs excluded in update.allele.info"
    stop(msg)
  }
  allele.info <- allele.info[id, ]
  allele.info
  
}
