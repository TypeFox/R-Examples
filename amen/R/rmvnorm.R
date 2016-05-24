#' Simulation from a multivariate normal distribution
#' 
#' Simulates a matrix where the rows are i.i.d. samples from a multivariate
#' normal distribution
#' 
#' 
#' @usage rmvnorm(n, mu, Sigma, Sigma.chol = chol(Sigma))
#' @param n sample size
#' @param mu multivariate mean vector
#' @param Sigma covariance matrix
#' @param Sigma.chol Cholesky factorization of \code{Sigma}
#' @return a matrix with \code{n} rows
#' @author Peter Hoff
#' @export rmvnorm
rmvnorm <-
function(n,mu,Sigma,Sigma.chol=chol(Sigma))
{
  # sample from a matrix normal distribution
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%Sigma.chol) +c(mu))
}
