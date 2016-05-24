#' Function to solve linear system using backsubsitution
#' 
#' Function to solve linear system using backsubsitution using Upper Triangular Matrix (Ux=c)
#' @param U     : Upper triangular matrix
#' @param c     : response
#' @return beta : Estiamted value 
optR.backsubsitution<-function(U, c) {  
  nROW<-nrow(U)
  beta<-matrix(rep(0, each=nROW), nrow=nROW, byrow=T)
  for(i in seq(nROW, 1, by=-1)) {
    if(U[i,i]!=0) {
      beta[i]<-(c[i] - apply(beta*U[i,], 2, sum))/U[i,i]
    } else
    {
      beta[i]<-0
    } 
  }
  return(beta)
}


#' Function to solve linear system using forward substitution
#' 
#' Function to solve linear system using backsubsitution using Upper Triangular Matrix (Ux=c)
#' @param L   : Lower triangular matrix
#' @param b   : response
#' @return y  : Estiamted value 
forwardsubsitution.optR<-function(L, b) {  
  nROW<-nrow(L)
  # Forward subsitution
  y<-b
  for(k in 2:nROW) {
    y[k]<-y[k]-sum(L[k,1:(k-1)]*y[1:(k-1)])
  }
  return(y)
}