globalARTP <-
function(data, B, K, gene_list, Gene="all", addit=FALSE,  covariable=NULL, family=binomial) {
  output <- NULL
  output$nPerm <- B
  if(Gene!="all") {
    colnames(gene_list) <- c("Id", "Gene")
    data <- data.frame(data$y, Selected_genes(Gene, gene_list, data))
    output$Gene <- Gene
  }
  Geno <- data[,-1]
  Trait <- factor(data[,1])
  if(K==1 || ncol(Geno) == 1) {
    output$Trunkpoint <- 1
    output$Kopt <- 1
    output$genevalue <- runPvalues(data, addit, covariable, family) 	
    return(output)
  }else{
    if(K>ncol(Geno)) K <- ncol(Geno)
    output$Trunkpoint <- K
    pvalors <- GeneratePvalues(data, B, addit, covariable, family)
    sb <- EstimatePvalue(B,K,pvalors)
    mP <- apply(sb,1,min)
    output$Kopt <- order(sb[1,])[1]
    output$genevalue <- sum(mP <= mP[1])/(B+1)
  }
  return(output)
}
