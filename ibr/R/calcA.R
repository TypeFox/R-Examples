calcA <- function(X,bx,kernelx="g") {
  X <- as.matrix(X)
  n <- nrow(X)
  KK <- matrix(0,ncol=n,nrow=n)
  Ddemi <- rep(0,n)
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
  for(i in 1:n){
    w <- poids(kernelx,X,bx,X[i,])
    KK[i,] <- w
    Ddemi[i] <- sqrt(1/sum(w))
  }
  A <- t(Ddemi*KK)*Ddemi
  return(list(A=A,Ddemi=Ddemi,df=sum(diag(A))))
}

