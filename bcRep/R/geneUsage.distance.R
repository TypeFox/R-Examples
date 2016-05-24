# Julia Bischof
# 12-08-2015

#library(vegan)
#library(proxy)


geneUsage.distance<-function(geneUsage.tab=NULL, names=NULL, method=c("bc","jaccard", "cosine"), cutoff=0){
  if(length(geneUsage.tab)==0){
    stop("--> gene usage table is missing")
  }
  if(length(method)!=1){
    stop("--> Please choose one method")
  }
  
  
  if(method=="bc"){
    dist<-"bray"
  }else if(method=="jaccard" && length(cutoff)==1){
    temp<-geneUsage.tab
    geneUsage.tab<-matrix(0, ncol=ncol(temp), nrow=nrow(temp))
    for(i in 1:nrow(geneUsage.tab)){
      geneUsage.tab[i,which(temp[i,]>cutoff)]<-1
    }
    dist<-"jaccard"
  }else{
    stop("--> Method is missing or unknown")
  }
  
  if(method=="cosine"){
    tab<-as.matrix(dist(geneUsage.tab, method="cosine"))
  }else{
    tab<-as.matrix(vegdist(geneUsage.tab, method = dist))
  }
  colnames(tab)<-if(length(names)>0){names}else{paste("Sample",1:nrow(geneUsage.tab),sep="")}
  rownames(tab)<-if(length(names)>0){names}else{paste("Sample",1:nrow(geneUsage.tab),sep="")}
  return(tab)  
}
