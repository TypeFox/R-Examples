#' Find Scale Parameter for Inverse Gamma Hyperprior 
#' 
#' This function implements a optimisation routine that computes the scale parameter \code{b}
#' of the inverse gamma prior for \eqn{\tau^2} when \eqn{a=b=\epsilon} with \eqn{\epsilon} small
#' for a given design matrix and prior precision matrix
#' such that approximately \eqn{P(|f(x_{k}|\le c,k=1,\ldots,p)\ge 1-\alpha}
#' When \code{a} unequal to \code{a} the shape parameter \code{a} has to be specified.



#' 
#' @param alpha denotes the 1-\eqn{\alpha} level.  
#' @param method with \code{integrate} as default.
#'        Currently no further method implemented.
#' @param Z the design matrix. 
#' @param c denotes the expected range of the function. 
#' @param eps denotes the error tolerance of the result, default is \code{.Machine$double.eps}. 
#' @param Kinv the generalised inverse of K. 
#' @param equals saying whether \code{a}=\code{b}. The default is FALSE due to the fact that a is a shape parameter.
#' @param a is the shape parameter of the inverse gamma distribution, default is 1.
#' @param type is either numerical integration (\code{integrate}) of to obtain the marginal distribution of \eqn{z_p'\beta}
#' 		  or the theoretical marginal t-distribution (\code{marginalt}). \code{marginalt} is the default value. 
#' @return an object of class \code{list} with values from \code{\link{uniroot}}.
#' @author Nadja Klein
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#' 
#' Stefan Lang and Andreas Brezger (2004). Bayesian P-Splines. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{13}, 183-212. 
#' @details Currently, the implementation only works properly for the cases \code{a} unequal \code{b}.
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
#' theta <- get_theta_ig(alpha = 0.01, method = "integrate", Z = Z, 
#'                       c = 3, eps = .Machine$double.eps, Kinv = Kinv, 
#'						 equals = FALSE, a = 1, type="marginalt")$root
#'


get_theta_ig <- function(alpha = 0.01, method = "integrate", Z, c = 3, eps = .Machine$double.eps, Kinv, equals = FALSE, a = 1,type="marginalt") 
  {
  warning("implementation currently does not work for equals=TRUE (a=b), when both values are small")
  if(method != "integrate")
    stop("method not existing")
  if(!(type %in% c("integrate","marginalt")))
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
  marginal_df <- function(f, lambda, ztz, a)
    {
	if(equals == TRUE)
	    a <- lambda
	if(type == "marginalt")  
	  {
	#  func <- function(b)
     #   {
	  df <- 2 * a
	  mu <- 0
	  sigma <- sqrt(ztz*(lambda/a))
      res <- pt((f-mu) / sigma, df = df) - pt((-f-mu) / sigma, df = df)  #- 1+ alpha
     #   }
	  #res <- try(uniroot(func, interval = c(eps4, 100)), TRUE)
	  #while(inherits(res, "try-error")) 
      #  {
      #  eps4 <- eps4 * 10
      #  res <- pt((c-mu) / sigma, df = df) - pt((-c-mu) / sigma, df = df)  - 1+ alpha#try(uniroot(func, interval = c(eps4, 100), Cov = ztKz, alpha = alpha), TRUE)
	  #  } 
	  } else {
	  integrand <- function(tau2, a=a) 
	    {
        dnorm(f, mean = 0, sd = sqrt(tau2 * ztz)) * ((lambda^a)/gamma(a) * tau2^(-a-1) * exp(-lambda/tau2))
        }
      res <- try(integrate(integrand, 0, Inf,a=a)$value, TRUE)
      while(inherits(res, "try-error")) 
	    {
        eps2 <- eps2 * 10
        res <- try(integrate(integrand, eps2, Inf,a=a)$value, TRUE)
        }
	  }
    return(res)
    }
  
  
  marginal_Pf <- function(lambda, Cov, alpha)
    {
    if(method == "integrate") 
	  {
      tempvar <- 0
      for(countnknots in 1:NROW(Z))
	    {
	    if(type=="integrate")
		  {
          contri <- try(2*integrate(Vectorize(marginal_df), -c, 0, lambda = lambda, ztz = Cov[countnknots], a=a)$value, TRUE)
          while(inherits(contri, "try-error")) 
		    {
            eps3 <- eps3 * 10
            contri <- try(2*integrate(Vectorize(marginal_df), -c, eps3, lambda = lambda, ztz = Cov[countnknots],a=a)$value, TRUE)
			}
          } else {
		  contri <- marginal_df(f=c, lambda = lambda, ztz = Cov[countnknots],a=a)
		  }
        tempvar <- tempvar + contri 
        }
      NROW(Z) - alpha - tempvar
      } else {
        stop("selected method not implemented.")
      }
    }

  result <- try(uniroot(marginal_Pf, interval = c(1000000000000*.Machine$double.eps, 1000), Cov = ztKz, alpha = alpha), TRUE)
   
  while(inherits(result, "try-error")) 
    {
    eps4 <- eps4 * 10
    result <- try(uniroot(marginal_Pf, interval = c(eps4, 1000), Cov = ztKz, alpha = alpha), TRUE)

    } 
  
  return(result)
}
