Wk<-function(x,tx,Pik,ck,b0=FALSE){
  
  if (b0 == TRUE){
    x<-as.matrix(cbind(1,x))}
  if (b0 == FALSE){
    x<-as.matrix(x)}
  
  tx<-as.matrix(tx)
  txpi<-as.matrix(t(x)%*%(1/Pik))
  V<-1/(Pik*ck)
  
  Wk<-(1/Pik)+((V*x)%*%solve(t(V*x)%*%x)%*%(tx-txpi))
  return(Wk)
}
