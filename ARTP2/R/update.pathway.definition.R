
update.pathway.definition <- function(pathway, exc.snps, exc.genes = NULL){
  
  if(length(exc.snps) == 0 && length(exc.genes) == 0){
    return(pathway)
  }
  
  id <- which(!(pathway$Gene %in% exc.genes))
  if(length(id) == 0){
    msg <- "All SNPs excluded in update.pathway.definition"
    stop(msg)
  }
  pathway <- pathway[id, ]
  
  id <- which(!(pathway$SNP %in% exc.snps))
  if(length(id) == 0){
    msg <- "All SNPs excluded in update.pathway.definition"
    stop(msg)
  }
  pathway <- pathway[id, ]
  rownames(pathway) <- NULL
  
  pathway
  
}
