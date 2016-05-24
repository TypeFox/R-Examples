"genotype.list"<-function(G, marker.type="MSW"){
gens<-list()
if(marker.type=="MSC" | marker.type=="SNP" | marker.type=="MSW"){
  if(length(G[1,])%%2!=0){stop("Genotypes have odd number of columns")}
  for(i in 1:(length(G[1,])/2)){
    gens[[i]]<-genotype(as.matrix(G[,((i*2)-1):(i*2)]))
    names(gens)[i]<-names(G[i*2]) 
  }
}
if(marker.type=="AFLP"){
  for(i in 1:length(G[1,])){
    gens[[i]]<-genotypeD(as.matrix(G[,i]))
    names(gens)[i]<-names(G[i]) 
  }
}
gens
}
