simulData <- function
### Function to simulate data \eqn{Y = X \beta + \sigma N(0, 1)}
(p=100, ##<< integer : number of variates. Should be >15 if \code{beta=NULL}
 n=100, ##<< integer : number of observations
 beta=NULL, ##<< vector with \code{p} components. See details.
 C=NULL, ##<<  matrix \code{p x p}. Covariance matrix of X. See details.
 r=0.95, ##<< scalar for calculating the covariance of X when \code{C=NULL}.
 rSN=10 ##<< scalar : ratio signal/noise
 ) {
  ##details<< When \code{beta} is \code{NULL}, then \code{p} should be
  ##greater than 15 and
  ##\code{beta=c(rep(2.5,5),rep(1.5,5),rep(0.5,5),rep(0,p-15))}
  if (is.null(beta)) {
    if (p<15) stop("Argument p should be greater than 15")
    beta <- c(rep(2.5,5),rep(1.5,5),rep(0.5,5),rep(0,p-15))
  }
  ##details<< When \code{C} is \code{NULL}, then \code{C} is block
  ##diagonal with \cr
  ## \code{C[a,b] = r**abs(a-b)} for \eqn{1 \le a, b \le 15} \cr
  ## \code{C[a,b] = r**abs(a-b)} for \eqn{16 \le a, b \le p} \cr
  if (is.null(C)) {
    C <- matrix(0,nrow=p,ncol=p)
    for (a in 1:15) {
      for (b in 1:15) {
        C[a,b] = r**abs(a-b)
      }
    }
    for (a in 16:p) {
      for (b in 16:p) {
        C[a,b] = r**abs(a-b)
      }
    }
  }
 # library(mvtnorm)
  ##note<< Library \code{mvtnorm} is loaded.
  ##details<< The lines of \code{X} are \code{n} i.i.d. gaussian variables with
  ##mean 0 and covariance matrix \code{C}.
  X <- rmvnorm(n=n,mean=rep(0,p),sigma=C)
  ##
  ##details<< The variance \code{sigma**2} equals the squared euclidean
  ##norm of \eqn{X \beta} divided by \code{rSN*n}.
  sigma <- sqrt(sum((X%*%beta)**2)/(rSN*n))
  Y <- X%*%beta + sigma*rnorm(n)
  return(
  ##value<< A list with components :
  list(
       Y=Y, ##<< vector \code{n} : \eqn{Y = X \beta + \sigma N(0, 1)}
       X=X,  ##<< matrix \code{n x p} : values of the covariates. See
       ##details.
       C=C,     ##<< matrix \code{p x p}. See details
       sigma=sigma, ##<< scalar. See details.
       beta=beta)) ##<< vector with \code{p} components. See details.
  ##end<<
} # fin simulData
# ---------------------------------------------------
