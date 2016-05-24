
update.deleted.snps <- function(deleted.snps, exc.snps, reason, comment){
  
  if(length(exc.snps) == 0){
    return(deleted.snps)
  }
  
  del.snps <- data.frame(SNP = exc.snps, reason = reason, comment = comment, stringsAsFactors = FALSE)
  deleted.snps <- rbind(deleted.snps, del.snps)
  dup <- duplicated(deleted.snps$SNP)
  deleted.snps <- deleted.snps[!dup, , drop = FALSE]
  
  deleted.snps
  
}
