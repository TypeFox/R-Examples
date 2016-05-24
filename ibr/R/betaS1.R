betaS1 <- function(n,U,tUy,eigenvaluesS1,ddlmini,k,lambda,Sgu,Qgu,index0) {
  prov <- rev(sumvalpr(k,n,rev(1-eigenvaluesS1),n-index0+1,n-ddlmini+1))
  beta <- U%*%(prov*tUy)
  qrSgu <- qr(Sgu)
  F2 <- qr.Q(qrSgu,complete=TRUE)[,-(1:ncol(Sgu))]
  ainv <- t(F2)%*%Qgu%*%F2
  diag(ainv) <- diag(ainv)+lambda
  Sp <- -lambda*F2%*%(solve(ainv))%*%t(F2)
  cgubeta <- as.vector(Sp%*%beta)/(-lambda)
  dgubeta <- solve(qr.R(qrSgu))%*%(t(qr.Q(qrSgu))%*%(as.matrix(beta)-Qgu%*%cgubeta))
  return(list(dgub=dgubeta,cgub=cgubeta))
}
