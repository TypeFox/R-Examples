#' Marginal Density for Given Scale Parameter and Half-Cauchy Prior for \eqn{\tau}
#' 
#' This function computes the marginal density of \eqn{z_p'\beta} for generalised beta prior hyperprior 
#' for \eqn{\tau^2} (half-Chauchy for \eqn{\tau})



#' 
#' @param f point the marginal density to be evaluated at. 
#' @param theta denotes the scale parameter of the generalised beta prior hyperprior for \eqn{\tau^2} (half-Chauchy for \eqn{\tau}).  
#' @param Z the row of the design matrix evaluated. 
#' @param Kinv the generalised inverse of K. 
#' @return the marginal density evaluated at point x.
#' @author Nadja Klein
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#' 
#' @import splines
#' @import GB2
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
#' mdf <- mdf_gbp(fgrid,theta=0.0028,Z=Z,Kinv=Kinv)
#'


mdf_gbp <- function(f, theta, Z, Kinv) 
  {
  ztKz <- diag(Z%*%Kinv%*%t(Z))
  eps2 <- .Machine$double.eps
  res <- c()
  for(countf in 1:length(f))
    {
    integrand <- function(tau2) 
	  {
      dnorm(f[countf], mean = 0, sd = sqrt(tau2 * ztKz)) * dgb2(tau2, shape1=1, scale=theta^2, shape2=0.5, shape3=0.5)
      }
    rescountf <- try(integrate(integrand, 0, Inf)$value, TRUE)
    while(inherits(rescountf, "try-error")) 
      {
      eps2 <- eps2 * 10
	  rescountf <- try(integrate(integrand, eps2, Inf)$value, TRUE)
	  }
	res <- c(res,rescountf)  
	}
  return(res)
  }
