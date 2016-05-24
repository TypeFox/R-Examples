#' Find Scale Parameter for Hyperprior
#' 
#' This function implements a optimisation routine that computes the scale parameter \eqn{\theta}
#' of the scale dependent hyperprior for a given design matrix and prior precision matrix
#' such that approximately \eqn{P(|f(x_{k}|\le c,k=1,\ldots,p)\ge 1-\alpha}



#' 
#' @param alpha denotes the 1-\eqn{\alpha} level.  
#' @param method either \code{integrate} or \code{trapezoid} with \code{integrate} as default.
#'        \code{trapezoid} is a self-implemented version of the trapezoid rule.
#' @param Z the design matrix. 
#' @param c denotes the expected range of the function. 
#' @param eps denotes the error tolerance of the result, default is \code{.Machine$double.eps}. 
#' @param Kinv the generalised inverse of K. 
#' @return an object of class \code{list} with values from \code{\link{uniroot}}.
#' @author Nadja Klein
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#' @import caTools
#' @import splines
#' @export
#' @examples
#' \dontrun{
#' 
#' set.seed(91179)
#' library(BayesX)
#' library(MASS)
#' # prior precision matrix to zambia data set
#' K <- read.gra(system.file("examples/zambia.gra", package="sdPrior"))
#' # generalised inverse of K
#' Kinv <- ginv(K)
#' 
#' # read data
#' dat <- read.table(system.file("examples/zambia_height92.raw", package="sdPrior"), header = TRUE)
#' 
#' # design matrix for spatial component
#' Z <- t(sapply(dat$district, FUN=function(x){1*(x==rownames(K))}))
#' 
#' # get scale parameter
#' theta <- get_theta(alpha = 0.01, method = "integrate", Z = Z, 
#'                             c = 3, eps = .Machine$double.eps, Kinv = Kinv)$root
#' } 


get_theta <- function(alpha = 0.01, method = "integrate", Z, c = 3, eps = .Machine$double.eps, Kinv) 
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
      dnorm(f, mean = 0, sd = sqrt(tau2 * ztz)) * dweibull(tau2, shape = 0.5, scale = lambda)
      }
    res <- try(integrate(integrand, 0, Inf)$value, TRUE)
    while(inherits(res, "try-error")) 
	  {
      eps2 <- eps2 * 10
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
      } else if(method == "sum") {
        stop("selected method not implemented yet.")
      } else if(method == "trapezoid") {
        tempvar <- 0
        for(countnknots in 1:NROW(Z)) 
		  {
          xseq <- 2 * (0:1000) * c / 1000 - c + 0.00001
          mdf <- sapply(xseq, FUN = marginal_df, lambda = lambda, ztz = Cov[countnknots])
          tempvar <- tempvar + alphafx[countnknots] + trapz(xseq, mdf)
          }
        1 - tempvar
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
