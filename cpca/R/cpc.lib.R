#-------------------------
# Main function `cpc`
#-------------------------

#' Function cpc. 
#'
#' This function computes the CPCA from a given set of covariance matrices 
#' (of different groups). 
#'
#' Currently, the only the stepwise algorithm by Trendafilov is supported.
#'
#' @name cpc
#' @param X An array of three dimensions: the 3rd dimension encodes the groups
#'   and the first two dimension contain the covariance matrices.
#' @param method The name of the method for computing the CPCA.
#'   The default value is \code{"stepwise"}, which is the stepwise algorithm by Trendafilov.
#' @param k The number of components to be computed (all if it is \code{0}).
#'   This paramter is valid if the given method supports 
#'   built-in ordering of the eigvenvectors.
#'   The default value is \code{0}, that means computing of all the components.
#' @param iter The maximum number of iterations.
#'   The parameter is valid for the stepwise algorithm by Trendafilov,
#'   that is applied in the power algorithm for estimation a single component.
#'   The default value is 30.
#' @param threshold The threshold value of the captured variance,
#'   which is reserved for further extensions.
#' @param ... Other parameters.
#' @return A list several slots: \code{CPC} rotation matrix with eigenvectors in columns;
#'   \code{ncomp} the number of components evaluated (equal to the number of columns in \code{CPC}).
#' @note This function adpats the original code in matlab written by Dr N. T. Trendafilov.
#' @references Trendafilov (2010). Stepwise estimation of common principal components. 
#'   Computational Statistics & Data Analysis, 54(12), 3446-3457. 
#'   doi:10.1016/j.csda.2010.03.010
#' @example inst/examples/function-cpc.R
#' @export
cpc <- function(X, method = "stepwise", k = 0, iter = 30, threshold = 0, ...)
{
  ### processing input argumets
  stopifnot(length(dim(X)) == 3)
  
  method <- match.arg(method, c("stepwise"))
  
  switch(method,
    stepwise = cpc_stepwise(X, apply(X, 3, nrow), k, iter, ...),
    stop("Error in swotch."))
}

cpc_stepwise <- function(X, n_g, k = 0, iter = 30, ...)
{
  p <- dim(X)[1]
  mcas <- dim(X)[3]

  # If k = 0 retrieve all components
  if(k == 0) {
    k <- p
  }
  
  # parameters: number of interations
  iter <- 15
  n <- n_g / sum(n_g)
  
  # output variables
  D <- array(0, dim=c(p, mcas))
  CPC <- array(0, dim=c(p, p))
  Qw <- diag(1, p)
  
  # components `s`
  s <- array(0, dim = c(p, p))
  for(m in 1:mcas) {
    s <- s + n[m]*X[, , m]
  }
  
  # variables
  res <- eigen(s)
  q0 <- res$vectors
  d0 <- diag(res$values, p)
  if(d0[1, 1] < d0[p, p]) {
    q0 <- q0[, ncol(q0):1]
  }
  
  # loop 'for ncomp=1:p'
  # Replaced by k so that only the first k components are retrieved
  for(ncomp in 1:k) {
   q <- q0[, ncomp]
   d <- array(0, dim=c(1, mcas))
   for(m in 1:mcas) {
    d[, m] <- t(q) %*% X[, , m] %*% q
   }
   
   # loop 'for i=1:iter'
   for(i in 1:iter) {
    s <- array(0, dim=c(p, p))
    for(m in 1:mcas) {
      s <- s + n_g[m] * X[, , m] / d[, m]
    }
           
    w <- s %*% q
    if( ncomp != 1) {
      w <- Qw %*% w
    }

    q <- w / as.numeric(sqrt((t(w) %*% w)))
    for(m in 1:mcas) {
      d[, m]  <- t(q) %*% X[, , m] %*% q
    }
    
   }
   # end of loop 'for i=1:iter'
   
   D[ncomp, ] <- d
   CPC[, ncomp] <- q
   Qw <- Qw - q %*% t(q)
  }
  # end of loop 'for ncomp=1:k'
  
  ### return
  out <- list(D = D[1:ncomp, ], CPC = CPC[, 1:ncomp], ncomp = ncomp)
  return(out)
}
