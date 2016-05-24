##' Fast pattern
##'
##' @title Fast pattern
##' @param x Matrix (binary) of patterns. Optionally if \code{y} is
##' also passed as argument, then the pattern matrix is defined as the
##' elements agreeing in the two matrices.
##' @param y Optional matrix argument with same dimensions as \code{x} (see above)
##' @param categories Default 2 (binary)
##' @param ... Optional additional arguments
##' @author Klaus K. Holst
##' @export
##' @examples
##' X <- matrix(rbinom(100,1,0.5),ncol=4)
##' fast.pattern(X)
##'
##' X <- matrix(rbinom(100,3,0.5),ncol=4)
##' fast.pattern(X,categories=4)
fast.pattern <- function(x,y,categories=2,...) {
  if (missing(y)) y <- NULL
   .Call("FastPattern",as.matrix(x),y,as.integer(categories))
}
 
