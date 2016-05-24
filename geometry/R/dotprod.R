##' If \code{x} and \code{y} are matrices, calculate the dot-product
##' along the first non-singleton dimension.  If the optional argument
##' \code{d} is given, calculate the dot-product along this
##' dimension. 
##' 
##' @title Compute the dot product of two vectors
##' @param x Matrix of vectors
##' @param y Matrix of vectors
##' @param d Dimension along which to calculate the dot product
##' @return Vector with length of \code{d}th dimension
##' @author David Sterratt
##' @keywords arith math array
##' @export
dot <- function(x, y, d=NULL) {
  if (is.vector(x)) 
    x <- matrix(x, ncol = 1)
  if (is.vector(y)) 
    y <- matrix(y, ncol = 1)
  ndim <- length(dim(x))
  
  ## Determine dimension along which to sum
  if (is.null(d)) {
    di <- which(dim(x) > 1)
    if (length(di == 0)) {
      d <- 1
    } else {
      d <- di[1]
    }
  }

  ## Check size of d
  if (d > ndim) {
    stop("d is larger than the number of dimensions of the data")
  }
  
  ## Use rowSums and colSums as they are more efficient than apply
  if (ndim == 2) {
    if (d == 2) {
      return(rowSums(x*y))
    }
    if (d == 1) {
      return(colSums(x*y))
    }
  }
  return(apply(x*y, d, sum))
}
