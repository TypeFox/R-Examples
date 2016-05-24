GSA.correlate=function(GSA.genesets.obj, genenames){
  nsets=length(GSA.genesets.obj$genesets)
  ngenes=unlist(lapply(GSA.genesets.obj$genesets,length))
  allgenes=unlist(GSA.genesets.obj$genesets)
  sets.in.exp=match(unique(allgenes),genenames)
  exp.in.sets=match(genenames,allgenes)
  
  cat("",fill=T)
  cat(c("Number of gene-sets:", nsets),fill=T)
  cat("",fill=T)
  cat(c("Total number of genes in gene-set collection:",sum(ngenes)),fill=T)
  cat(c("Total number of unique genes in gene-set collection:",length(unique(allgenes))),fill=T)
  cat("",fill=T)
   cat(c("Total number of genes in  genenames list:",length(genenames)),fill=T)
  cat(c("Total number of unique genes in genenames list:",length(unique(genenames))),fill=T)
  cat("",fill=T)
cat(c("Number of unique genes in both collections:",sum(!is.na(sets.in.exp))),fill=T)

  nn=rep(NA,nsets)
  for(i in 1:nsets){
    nn[i]=sum(!is.na(match(GSA.genesets.obj$genesets[[i]],genenames)))
  }
  cat("",fill=T)
  cat(c("Quantiles of fraction coverage of gene-sets"),fill=T)
  print(quantile(nn/ngenes, seq(0,1,by=.1)),digits=4)
  return()
}
