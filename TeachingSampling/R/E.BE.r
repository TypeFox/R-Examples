E.BE<-function(y,prob){
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "N"
  Total<-matrix(NA,nrow=4,ncol=dim(y)[2])
  rownames(Total)=c("Estimation", "Standard Error","CVE","DEFF")
  colnames(Total)<-names(y)
  
  for(k in 1:dim(y)[2]){
    ty<-sum(y[,k])/prob
    Vty<-(1/prob)*((1/prob)-1)*sum(y[,k]^2)
    CVe<-100*sqrt(Vty)/ty
    n<-length(y[,k])
    N<-n/prob
    VMAS<-(N^2)*(1-(n/N))*var(y[,k])/(n)
    DEFF<-Vty/VMAS
    Total[,k]<-c(ty,sqrt(Vty),CVe,DEFF)
  }
  return(Total)
}