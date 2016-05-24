poids <- function(kernelx,X,bx,valx,n,p){
  VALX <- matrix(rep(valx,n),ncol=p,byrow=TRUE)
  BX <- matrix(rep(bx,n),ncol=p,byrow=TRUE)
    MAT <- (X-VALX)/BX
  if(kernelx=="e"){noyau <- epane}
  if(kernelx=="g"){noyau <- gaussien}
  if(kernelx=="q"){noyau <- quartic}
  if(kernelx=="u"){noyau <- uniform}
  vv <- noyau(MAT)
  vv <- apply(vv,1,prod)
  return(W=vv)
}
tracekernel <- function(X,bx,kernelx,n,p){
  H <- matrix(0,ncol=n,nrow=n)
    for(i in 1:n){
      w <- poids(kernelx,X,bx,X[i,],n,p)
      H[i,] <- w/sum(w)
    }
  trace <- sum(diag(H))
  return(trace=trace)
}
departnoyau <- function(df,x,kernel,dftobwitmax,n,p,dfobjectif) {
  bandwidth <- bwchoice(x,df,kernel,dftobwitmax)
  resultat <- tracekernel(x,bandwidth,kernel,n,p)-dfobjectif
  return(resultat)
}
