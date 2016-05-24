#' @title Test if an object has one-dimension
#' 
#' @description
#' Returns \code{TRUE} if an object is a vector or a one-dimensional
#' matrix, \code{FALSE} otherwise
#' 
#' @param x an R object
#' @return whether x is one-dimensional
#' @export
#' @examples
#' # vector
#' is_one_dim(1:5)  # TRUE
#' 
#' # factor
#' is_one_dim(iris$Species)  # TRUE
#' 
#' # one row matrix
#' is_one_dim(matrix(1:5, 1, 5))  # TRUE
#' 
#' # one column matrix
#' is_one_dim(matrix(1:5, 5, 1))  # TRUE
#' 
#' # general matrix (nrow>1, ncol>1)
#' is_one_dim(matrix(1:9, 3, 3))  # FALSE
#' 
#' # general data frame
#' is_one_dim(iris)  # FALSE
is_one_dim <- function(x)
{
  one_dim = TRUE
  if (lacks_dim(x)) {
    if (is.list(x)) one_dim = FALSE
  } else {
    if (dim(x)[1L] > 1 && dim(x)[2L] > 1)
      one_dim = FALSE
  }
  # output
  one_dim
}
