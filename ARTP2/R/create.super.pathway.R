
create.super.pathway <- function(pathway){
  
  super.pathway <- NULL
  npath <- length(pathway)
  header <- c("SNP", "Gene", "Chr")
  
  for(i in 1:npath){
    super.pathway <- rbind(super.pathway, pathway[[i]][, header])
  }
  
  super.pathway <- super.pathway[!duplicated(super.pathway), ]
  super.pathway
  
}
