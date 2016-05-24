kernelSx <- function(kernelx="g",X,Xetoile,bx){
  poids <- function(kernelx="g",X,bx,valx){
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
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
  X <- as.matrix(X)
  n <- nrow(X)
  Xetoile <- as.matrix(Xetoile)
  netoile <- nrow(Xetoile)
  H <- matrix(0,ncol=n,nrow=netoile)
  for (i in 1:netoile){
    w <- poids(kernelx,X,bx,Xetoile[i,])
    H[i,] <- w/sum(w)
  }
 return(H)
}
