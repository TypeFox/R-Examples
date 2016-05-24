E.Beta<-function(N, n, y, x, ck=1, b0=FALSE){
  if (b0 == TRUE) {
    x<-as.data.frame(cbind(1,x))
  }
  #---------------
  Q <- dim(as.matrix(y))[2]
  P <- dim(as.matrix(x))[2]
  #---------------
  Total<-array(NA,c(3,P,Q))
  rownames(Total)=c("Beta estimation", "Standard Error","CVE")
  colnames(Total)<-names(x)
  dimnames(Total)[[3]]<-names(y)
  #---------------
  Pik <- rep(n/N,n)
  for(q in 1:Q){
    yq<-as.matrix(y[,q])
    x<-as.matrix(x)
    ck<-as.numeric(unlist(ck))
    V<-1/(Pik*ck)
    bq<-solve(t(V*x)%*%x)%*%(t(V*x)%*%yq)
    ek <- yq - x%*%bq
    uk <- c(ek)*x 
    Varuk <- (N^2/n)*(1-(n/N))*var(uk)
    P1 <- solve(t(V*x)%*%x)
    Vbeta <- as.matrix(P1)%*%as.matrix(Varuk)%*%as.matrix(P1)
    Vbeta <- diag(Vbeta)
    CVe <-100*sqrt(Vbeta)/bq  
    #---------------  
    if(Q == 1){
      Total[1,,]<-bq
      Total[2,,]<-sqrt(Vbeta)
      Total[3,,]<-CVe  
    }
    if(P == 1 & Q > 1){
      Total[1,,][q]<-bq
      Total[2,,][q]<-sqrt(Vbeta)
      Total[3,,][q]<-CVe
    }
    if(Q > 1 & P > 1){
      Total[1,,][,q]<-bq
      Total[2,,][,q]<-sqrt(Vbeta)
      Total[3,,][,q]<-CVe  
    }
    #---------------
  }
  return(Total)
}