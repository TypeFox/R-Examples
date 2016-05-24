"modeG"<-function(postG, threshold=0){
  id<-rownames(postG[[1]])
  G<-as.list(1:length(postG))
  names(G)<-names(postG)
  n<-sum(postG[[1]][1,])
  for(i in 1:length(G)){
    an<-allele.names(as.genotype(colnames(postG[[i]])))
    G[[i]]<-as.genotype(colnames(postG[[i]])[apply(postG[[i]],1,which.max)], alleles=an)
    G[[i]][which(apply(postG[[i]],1,max)/n < threshold)]<-NA
  }
  list(G=G, id=id)
}
 
