Cov.Gauss <- function(d.Beta,d.sigma,sigma,rs0,delta,X,cu,const)
{
  n   <- length(rs0)
  p   <- ncol(X)
  JAC <- matrix(0, nrow=(p+1), ncol=(p+1))
  J11 <- TMLjac11.G(d.Beta,d.sigma,rs0,delta,X,cu)
  J12 <- TMLjac12.G(d.Beta,d.sigma,rs0,delta,X,cu)
  J21 <- TMLjac21.G(d.Beta,d.sigma,rs0,delta,X,cu)
  J22 <- TMLjac22.G(d.Beta,d.sigma,rs0,delta,X,cu)
  JAC[1:p,1:p]     <- J11
  JAC[1:p,(p+1)]   <- J12
  JAC[(p+1),1:p]   <- J21
  JAC[(p+1),(p+1)] <- J22
  JAC <- JAC/sigma
  JIC <- solve(JAC)
  Q   <- QMatrix.G(d.Beta,d.sigma,rs0,delta,X,cu,const)
  Cov <- t(JIC) %*% Q %*% JIC /n
  Cov
}