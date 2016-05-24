getPvalBH<-function(x)
 {
  pvals<-as.numeric(unlist(x))
  p.adjust.M<-p.adjust.methods
  p.adj<-sapply(p.adjust.M, function(meth) p.adjust(pvals,meth))
  dd<-data.frame(pval=pvals,pval.BH=p.adj[,5])
  o<-order(dd[,2])
  ddEnd<-data.frame(region=o,dd[o,])
  dimnames(ddEnd)[[1]]<-c(1:nrow(ddEnd))
  ddEnd
 }
