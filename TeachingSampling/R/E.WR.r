E.WR<-function(N,m,y){
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "N"
  Total<-matrix(NA,nrow=4,ncol=dim(y)[2])
  rownames(Total)=c("Estimation", "Standard Error","CVE","DEFF")
  colnames(Total)<-names(y)
  
  for(k in 1:dim(y)[2]){
    
    ty<-(N/m)*sum(y[,k])
    Vty<-(N^2/m)*var(y[,k])
    CVe<-100*sqrt(Vty)/ty
    VMAS<-(N^2)*(1-(m/N))*var(y[,k])/(m)
    DEFF<-Vty/VMAS
    Total[,k]<-c(ty,sqrt(Vty),CVe,DEFF)
  }
  return(Total)
}
