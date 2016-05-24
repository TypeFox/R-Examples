T.SIC<-function(y,Cluster){
  
  Cluster<-as.factor(Cluster)
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "Ni"
  
  nI<-length(levels(Cluster))
  
  Total<-matrix(NA,nrow=nI,ncol=dim(y)[2],)
  rownames(Total)<-levels(Cluster)
  colnames(Total)<-names(y)
  Cluster<-as.factor(as.integer(Cluster))
  
  for(k in 1: nI){
    e<-which(Cluster==k)
    ye<-y[e,]
    ye<-as.matrix(ye)
    tye<-colSums(ye)
    Total[k,]<-tye
  }
  Total<-as.matrix(Total)
  return(Total)
}