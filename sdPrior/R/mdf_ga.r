#' Marginal Density for Given Scale Parameter and Half-Normal Prior for \eqn{\tau}
#' 
#' This function computes the marginal density of \eqn{z_p'\beta} for gamma priors for \eqn{\tau^2}
#' (referring to a half-normal prior for \eqn{\tau}).



#' 
#' @param f point the marginal density to be evaluated at. 
#' @param theta denotes the scale parameter of the gamma hyperprior for \eqn{\tau^2} (half-normal for \eqn{\tau}).  
#' @param Z the row of the design matrix evaluated. 
#' @param Kinv the generalised inverse of K. 
#' @return the marginal density evaluated at point x.
#' @author Nadja Klein
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#' 
#' @import splines
#' @export
#' @examples
#' set.seed(123)
#' library(MASS)
#' # prior precision matrix (second order differences) 
#' # of a spline of degree l=3 and with m=20 inner knots
#' # yielding dim(K)=m+l-1=22
#' K <- t(diff(diag(22), differences=2))%*%diff(diag(22), differences=2)
#' # generalised inverse of K
#' Kinv <- ginv(K)
#' # covariate x
#' x <- runif(1)
#' Z <- matrix(DesignM(x)$Z_B,nrow=1)
#' fgrid <- seq(-3,3,length=1000)
#' mdf <- mdf_ga(fgrid,theta=0.0028,Z=Z,Kinv=Kinv)
#'


mdf_ga <- function(f, theta, Z, Kinv) 
  {
  ztKz <- diag(Z%*%Kinv%*%t(Z))
  cp <- 1/(theta*pi*sqrt(ztKz))
  arg <- sqrt(f^2/(theta^2*ztKz))
  res <- cp*besselK(x=arg,nu=0)
  return(res)
}
