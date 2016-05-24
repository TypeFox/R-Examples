summary.GSA.genesets=function(object,...){
  nsets=length(object$genesets)
  ngenes=unlist(lapply(object$genesets,length))
  allgenes=unlist(object$genesets)
  cat("",fill=T)
  cat(c("Number of gene-sets:", nsets),fill=T)
  cat("",fill=T)
  cat(c("Total number of genes in gene-set collection:",sum(ngenes)),fill=T)
  cat(c("Total number of unique genes in gene-set-collection:",length(unique(allgenes))),fill=T)
  cat("",fill=T) 
  cat(c("Table of numbers of genes in gene-sets:"),fill=T)
  tt=table(ngenes,dnn=NULL)
k=15
  nn=trunc(length(tt)/k)
 cat("",fill=T)
for(i in 1:(nn+1)){
  i1=1+k*(i-1)
  i2=min(i1+k-1, length(tt))
  print(tt[i1:i2],dnn=NULL)
  cat("",fill=T)
}
 
  return()
}
