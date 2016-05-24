globalSimes <-
function(data, B, gene_list, Gene="all", addit=FALSE, covariable=NULL, family=binomial) {
  Geno <- data[,-1]
  Trait <- factor(data[,1])
  output <- NULL
  output$nPerm <- B
  if(Gene!="all") {
    data <- data.frame(Trait, Selected_genes(Gene, gene_list, data))
    output$Gene <- Gene  
  }
  x <- GeneratePvalues(data, B, addit, covariable, family)
  output$genevalue <- min(length(x)*sort(x)/seq(1,length(x), by=1))
  return(output)
}
