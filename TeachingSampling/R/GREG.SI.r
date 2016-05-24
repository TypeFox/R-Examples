GREG.SI<-function(N,n,y,x,tx,b,b0=FALSE){
  y<-as.data.frame(y)
  x<-as.matrix(x)
  pik<-rep(n/N,n)
  dk<-1/pik
  if (b0 == TRUE){
    x<-as.matrix(cbind(1,x))}
  
  Total<-matrix(NA,nrow=3,ncol=dim(y)[2])
  rownames(Total)=c("Estimation", "Standard Error","CVE")
  colnames(Total)<-names(y)
  
  for(k in 1:dim(y)[2]){
    
    xHT <-  t(x)%*%dk
    yHT <- sum(y[,k]*dk)
    ty <- yHT + (tx-t(xHT))%*%as.matrix(b[,k])
    e <- y[,k]-(x%*%as.matrix(b[,k]))
    Vty <- (N^2)*(1-(n/N))*var(e)/(n)
    CVe <- 100*sqrt(Vty)/ty
    Total[,k] <- c(ty,sqrt(Vty),CVe)
  }
  return(Total)
}