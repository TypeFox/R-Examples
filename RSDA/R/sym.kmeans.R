sym.kmeans <-
function(sym.data,k=3,iter.max=10,nstart=1,
                     algorithm=c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen")) {
  algorithm<-match.arg(algorithm)
  idn<-all(sym.data$sym.var.types=='$I')
  if(idn==FALSE)
    stop("The two variables have to be interval type")         
    nn<-sym.data$N
    mm<-sym.data$M
    centers<-matrix(0,nn,mm)
    centers<-as.data.frame(centers)
    rownames(centers)<-sym.data$sym.obj.names
    colnames(centers)<-sym.data$sym.var.names
    for(i in 1:nn) 
      for(j in 1:mm)
        centers[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                         sym.var(sym.data,j)$var.data.vector[i,2])/2
    return(kmeans(centers,k,iter.max,nstart,algorithm))  
}
