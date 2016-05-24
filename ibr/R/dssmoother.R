dssmoother <- function(X,Y=NULL,lambda,m,s) {
  n <- nrow(X)
  d <- ncol(X)
  Sgu <-  DuchonS(X,m)
  qrSgu <- qr(Sgu)
  F2 <- qr.Q(qrSgu,complete=TRUE)[,-(1:ncol(Sgu))]
  Qgu <- DuchonQ(X,0,m,s,symmetric=TRUE)
  ainv <- t(F2)%*%Qgu%*%F2
  diag(ainv) <- diag(ainv)+lambda
  Sp <- -lambda*F2%*%(solve(ainv))%*%t(F2)
  if (! is.null(Y)) {
    cgu <- as.vector(Sp%*%Y)/(-lambda)
    dgu <- solve(qr.R(qrSgu))%*%(t(qr.Q(qrSgu))%*%(as.matrix(Y)-Qgu%*%cgu))
    diag(Sp) <- rep(1,n)+diag(Sp)
    return(list(H=Sp,Sgu=Sgu,Qgu=Qgu,dgu=dgu,cgu=cgu))
  }  else {
    diag(Sp) <- rep(1,n)+diag(Sp)
    return(list(H=Sp,Sgu=Sgu,Qgu=Qgu,dgu=NULL,cgu=NULL))
  }
}
