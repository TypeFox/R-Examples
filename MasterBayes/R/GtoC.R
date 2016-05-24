"GtoC"<-function(G, biallelic=FALSE){
  
  for(i in 1:length(G)){
    if(is.genotype(G[[i]]) & biallelic==FALSE){
     G[[i]]<-matrix(match(allele(G[[i]]),allele.names(G[[i]]), nomatch=-998), length(G[[i]]),2)-1
    }
    if(is.genotype(G[[i]]) & biallelic==TRUE){
     G[[i]]<-matrix(match(allele(G[[i]]),allele.names(G[[i]]), nomatch=NA), length(G[[i]]),2)-1
     G[[i]]<-as.matrix(G[[i]][,1]+G[[i]][,2])
     G[[i]][which(is.na(G[[i]])==TRUE)]<--999
    }
    if(is.genotypeD(G[[i]])){
     G[[i]]<-as.numeric(G[[i]])-1
     G[[i]][which(is.na(G[[i]])==TRUE)]<--999
    }
  }
  if(biallelic==FALSE){
    G<-c(t(matrix(unlist(G), length(G[[1]][,1]),2*length(G))))
  }else{
    G<-c(t(matrix(unlist(G), length(G[[1]]),length(G))))
  }
  G
}
