
update.deleted.genes <- function(deleted.genes, exc.genes, reason){
  
  if(length(exc.genes) == 0){
    return(deleted.genes)
  }
  
  del.genes <- data.frame(Gene = exc.genes, reason = reason, stringsAsFactors = FALSE)
  deleted.genes <- rbind(deleted.genes, del.genes)
  dup <- duplicated(deleted.genes$Gene)
  deleted.genes <- deleted.genes[!dup, , drop = FALSE]
  
  deleted.genes
  
}

