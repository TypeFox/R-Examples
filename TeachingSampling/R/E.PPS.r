E.PPS<-function(y,pk){
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "N"
  Total<-matrix(NA,nrow=4,ncol=dim(y)[2])
  rownames(Total)=c("Estimation", "Standard Error","CVE","DEFF")
  colnames(Total)<-names(y)
  m<-length(pk)
  for(k in 1:dim(y)[2]){
    ty<-sum(y[,k]/pk)/m
    Vty<-(1/m)*(1/(m-1))*sum((y[,k]/pk-ty)^2)
    CVe<-100*sqrt(Vty)/ty 
    N<-(1/m)*sum(1/pk)
    VMAS<-(N^2)*(1-(m/N))*var(y[,k])/(m)
    DEFF<-Vty/VMAS
    Total[,k]<-c(ty,sqrt(Vty),CVe,DEFF)
  }
  return(Total)
}