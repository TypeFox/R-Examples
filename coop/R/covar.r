#' Covariance
#' 
#' An optimized, efficient implemntation for computing covariance.
#' 
#' @details
#' See \code{?coop-package} for implementation details.
#' 
#' @param x
#' A numeric matrix or vector.
#' @param y
#' A vector (when \code{x} is a vector) or missing (blank) when 
#' \code{x} is a matrix.
#' @param use
#' The NA handler, as in R's \code{cov()} and \code{cor()}
#' functions.  Options are "everything", "all.obs", and 
#' "complete.obs".
#' 
#' @return
#' The covariance matrix.
#' 
#' @examples
#' x <- matrix(rnorm(10*3), 10, 3)
#' 
#' coop::pcor(x)
#' coop::pcor(x[, 1], x[, 2])
#' 
#' @author Drew Schmidt
#' @seealso \code{\link{cosine}}
#' @export
covar <- function(x, y, use="everything") UseMethod("covar")



#' @export
covar.matrix <- function(x, y, use="everything")
{
  co_matrix(x, y, CO_VAR, use)
}



#' @export
covar.default <- function(x, y, use="everything")
{
  co_vecvec(x, y, CO_VAR, use)
}
