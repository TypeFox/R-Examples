#R functions to deal with multivariate normals efficiently using Rcpp and RcppEigen functions

#draw n samples from a multivariate normal distribution with covariance matrix S and mean mu
rcpp_rmvnorm <- function(n,S,mu)
{
  chol_S <- rcppeigen_get_chol(S)
  m <- dim(chol_S)[1]
  mat <- matrix(rnorm(n*m),nrow=m,ncol=n)
  return (t(chol_S%*%mat+mu))
}
#draw n samples from a multivariate normal distribution with covariance matrix S and mean mu with a more stable version of Cholesky
rcpp_rmvnorm_stable <- function(n,S,mu)
{
  if(length(mu) == 1){
    return(rnorm(n,mean=mu,sd=sqrt(S)))
  }
  else{
    chol_S <- rcppeigen_get_chol_stable(S)
    diag_s <- rcppeigen_get_chol_diag(S)
    diag_s[which(diag_s < 0)] <- 0
    diag_s <- diag(sqrt(diag_s))
    m <- dim(chol_S)[1]
    mat <- matrix(rnorm(n*m),nrow=m,ncol=n)
    return (t(chol_S%*%diag_s%*%mat+mu))
  }
}
#invert a Toeplitz matrix (useful for GPs where points are evenly spaced)
tinv <- function(A){
  M <- A/(A[1,1])
  N <- dim(M)[1]
  r <- M[1,2:N]
  y <- durbin(r,N-1)
  return(trench(r,y,N)/A[1,1])
}
#evaluate the log density of a multivariate normal distribution with covariance S and mean mu at point x, with a flag for whether the covariance matrix is Toeplitz
rcpp_log_dmvnorm <- function(S,mu,x, istoep)
{
  n <- length(x)
  diag <- rcppeigen_get_diag(S)
  if (istoep == TRUE)
  {
    inv_S <- tinv(S)
  }
  else
  {
    inv_S <- rcppeigen_invert_matrix(S)
  }
  return((-n/2)*log(2*pi)-sum(log(diag))-(1/2)*(x-mu)%*%inv_S%*%(x-mu))
}
