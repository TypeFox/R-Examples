#' Decompose covariance matrix
#' 
#' \code{decomp.cov} decomposes a covariance matrix \code{v}.  If \code{A = decomp.cov(v)}, then \code{tcrossprod(A, A) == v}. 
#' 
#' The \code{"chol"} method is the fastest, but must unstable.  The \code{"eigen"} method is slower, but more stable.  The \code{"svd"} method is the slowest method, but should be the most stable. 
#' 
#' @param v An \eqn{n \times n} covariance matrix.
#' @param method The method used to decompose \code{v}. valid options are \code{"chol"}, \code{"eigen"}, or \code{"svd"}.
#' 
#' @return Returns an \eqn{n \times n} matrix.
#' 
#' @author Joshua French
#' @export
#' @examples 
#' # generate data
#' n = 100
#' coords = matrix(runif(n*2), nrow = n, ncol = 2)
#' d = as.matrix(dist(coords))
#' # create covariance matrix
#' v = 3*exp(-d/2) + 0.1*diag(n)
#' 
#' # decompose v using the three methods
#' d1 = decomp.cov(v, "chol")
#' d2 = decomp.cov(v, "eigen")
#' d3 = decomp.cov(v, "svd")
#' 
#' # verify accuracy of decompositions
#' range(v - tcrossprod(d1))
#' range(v - tcrossprod(d2))
#' range(v - tcrossprod(d3))

decomp.cov <- function(v, method = "eigen")
{
  #sanity check
  if(!(is.matrix(v)  || inherits(v, "Matrix")))
  {
    stop("v should be a matrix or Matrix")
  }
  if(nrow(v)!=ncol(v))
  { stop("v must be a square numeric matrix")}
  if(!is.element(method, c("chol", "eigen", "svd")))
  { stop("method must be 'chol', 'eigen', or 'svd'")}
  
  if(method == "eigen")
  {
    eigenv <- eigen(v)
    return(eigenv$vectors %*% diag(sqrt(pmax(eigenv$values,0))))
  }else if(method == "chol")
  { 
    return(t(chol(v))) 
  }else if(method == "svd")
  {
    svdv <- svd(v)
    return(tcrossprod(svdv$u %*% diag(sqrt(svdv$d)), svdv$v))
  }
}
