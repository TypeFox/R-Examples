

  #######################
  #### function GtoD2 ####
  #######################
  
  
 GtoD2 <- function(G){
  
  if (class(G)!="Gram")
   stop(" 'G' must be of class Gram")
   
  n<-ncol(G)
  onesn <- matrix(rep(1,n),nrow=n) #ones vector
  
  Delta<-t(t(diag(G)))%*%t(onesn)+onesn%*%t(diag(G))-2*G
  class(Delta)<-"D2"
  return(Delta)
 }