E.SI<-function(N,n,y){
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "N"
  Total<-matrix(NA,nrow=4,ncol=dim(y)[2])
  rownames(Total)=c("Estimation", "Standard Error","CVE","DEFF")
  colnames(Total)<-names(y)
  pik<-matrix(n/N, nrow=n, ncol=1)
  dk<-1/pik
  for(k in 1:dim(y)[2]){
    ty<-sum(y[,k]*dk)
    Vty<-(N^2)*(1-(n/N))*var(y[,k])/(n)
    CVe<-100*sqrt(Vty)/ty
    VMAS<-(N^2)*(1-(n/N))*var(y[,k])/(n)
    DEFF<-Vty/VMAS
    Total[,k]<-c(ty,sqrt(Vty),CVe,DEFF)
  }
  return(Total)
}