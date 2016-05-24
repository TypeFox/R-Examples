#' Find Scale Parameter for Hyperprior for Variances Where the Standard Deviations have an 
#' Approximated (Differentiably) Uniform Distribution. 
#' 
#' This function implements a optimisation routine that computes the scale parameter \eqn{\theta}
#' of the prior \eqn{\tau^2} (corresponding to a differentiably approximated version of the uniform prior for \eqn{\tau}) 
#' for a given design matrix and prior precision matrix
#' such that approximately \eqn{P(|f(x_{k}|\le c,k=1,\ldots,p)\ge 1-\alpha}



#' 
#' @param alpha denotes the 1-\eqn{\alpha} level.  
#' @param method with \code{integrate} as default.
#'        Currently no further method implemented.
#' @param Z the design matrix. 
#' @param c denotes the expected range of the function. 
#' @param eps denotes the error tolerance of the result, default is \code{.Machine$double.eps}. 
#' @param Kinv the generalised inverse of K. 
#' @return an object of class \code{list} with values from \code{\link{uniroot}}.
#' @author Nadja Klein
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#' 
#' Andrew Gelman (2006). Prior Distributions for Variance Parameters in Hierarchical Models. 
#' \emph{Bayesian Analysis}, \bold{1}(3), 515--533. 
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
#' theta <- get_theta_aunif(alpha = 0.01, method = "integrate", Z = Z, 
#'                             c = 3, eps = .Machine$double.eps, Kinv = Kinv)$root
#'


get_theta_aunif <- function(alpha = 0.01, method = "integrate", Z, c = 3, eps = .Machine$double.eps, Kinv) 
  {
  if(method != "integrate")
    stop("method not existing")
  ztKz <- diag(Z%*%Kinv%*%t(Z))
  #number of grids of x are given by the rows of Z
  if(NROW(Z) == 1 | NCOL(Z) == 1) 
    {
    nknots <- 1
    } else {
    nknots <- NROW(Z)
    }
  #weights such that sum(weights) = nknots
  weights <- rep(1, nknots)
  #alpha-level for each point x (alphax = weight * alpha / nknots)
  alphafx <- alpha * weights / nknots
  eps2 <- eps3 <- eps4 <- eps
  marginal_df <- function(f, lambda, ztz)
    {
    integrand <- function(tau2) 
	  {
      dnorm(f, mean = 0, sd = sqrt(tau2 * ztz)) * dapprox_unif(tau2,scale=lambda)#,tildec=6.931472)
      }
	lowval <- eps2
    res <- try(integrate(integrand, lowval, Inf)$value, TRUE)
    while(inherits(res, "try-error")) 
	  {
      eps2 <- lowval * (1.01)
      res <- try(integrate(integrand, eps2, Inf)$value, TRUE)
      }
    return(res)
    }
  
  
  marginal_Pf <- function(lambda, Cov, alpha)
    {
    if(method == "integrate") 
	  {
      tempvar <- 0
      for(countnknots in 1:NROW(Z)) {
        contri <- try(2*integrate(Vectorize(marginal_df), -c, 0, lambda = lambda, ztz = Cov[countnknots])$value, TRUE)
        while(inherits(contri, "try-error")) {
          eps3 <- eps3 * 10
          contri <- try(2*integrate(Vectorize(marginal_df), -c, eps3, lambda = lambda, ztz = Cov[countnknots])$value, TRUE)
        }
        tempvar <- tempvar + contri 
      }
      NROW(Z) - alpha - tempvar
      } else {
        stop("selected method not implemented.")
      }
    }
  result <- try(uniroot(marginal_Pf, interval = c(1000000000000*.Machine$double.eps, 1000), Cov = ztKz, alpha = alpha), TRUE)
  if(inherits(result, "try-error"))
    {
	if(grepl("not of opposite sign",as.character(result,"condition")))
	  result <- try(uniroot(marginal_Pf, interval = c(0.118, 1000), Cov = ztKz, alpha = alpha), TRUE)
	}
  if(inherits(result, "try-error"))
    {
	if(grepl("not of opposite sign",as.character(result,"condition")))
	  {
	  findlowercands <- seq(0.12,1,length=100)
      vals <- c()
	  for(counti in 1:length(findlowercands))
	    {
	    vals <- c(vals,(2*NROW(Z)*integrate(Vectorize(marginal_df), -c, 0, lambda = findlowercands[counti], ztz = ztKz[1])$value-(NROW(Z)-alpha)))
	    }
	    index <- which(vals>=0)	
	    result2 <- list(warns="finding the root is not possible", root=findlowercands[which.min(vals[index])])	
	    return(result2)
	  } else {
	  warning("numerical optimisation instable, scale parameter might be wrong and computation may take some time")
	  findlower <- 0.12
	  while(((2*NROW(Z)*integrate(Vectorize(marginal_df), -c, 0, lambda = findlower, ztz = ztKz[1])$value-(NROW(Z)-alpha))<=0) && findlower>=0.000001)
        {
	    findlower <- findlower*0.99
	    }
      result <- try(uniroot(marginal_Pf, interval = c(findlower, 1000), Cov = ztKz, alpha = alpha), TRUE)      
	  }
	}   
  return(result)
}
