#' Extract the terms of a multivariate polynomial.
#' 
#' Compute the terms of an mpoly object as a mpolyList.
#' 
#' @param x an object of class mpoly
#' @param ... additional parameters
#' @return An object of class mpolyList.
#' @usage \method{terms}{mpoly}(x, ...)
#' @export
#' @examples
#' 
#' \dontrun{ .Deprecated issues a warning
#' 
#' x <- mp("x^2 - y + x y z")
#' terms(x)
#' monomials(x)
#' 
#' }
#' 
terms.mpoly <- function(x, ...){
  .Deprecated("monomials")
  # mpolyList <- lapply(x, function(v){
  #   l <- list(v)
  #   class(l) <- 'mpoly'
  #   l
  # })
  # class(mpolyList) <- 'mpolyList'
  # mpolyList
  monomials(x)
}
