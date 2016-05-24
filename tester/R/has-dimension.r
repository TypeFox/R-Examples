#' @title Has dimension?
#' @description 
#' \code{has_dimension} and \code{has_dim} test if an object 
#' has dimension (i.e. \code{dim}) \cr
#' \code{lacks_dimension} and \code{lacks_dim} test if an object 
#' lacks dimension
#' 
#' @param x an R object
#' @aliases has_dimension has_dim lacks_dimension lacks_dim
#' @export has_dimension has_dim lacks_dimension lacks_dim
#' @examples
#' m = matrix(1:12, 4, 3)
#' a = as.array(letters)
#' has_dim(m) # TRUE
#' has_dimension(a)
#' 
#' has_dimension(iris) # TRUE
#' 
#' has_dim(matrix(1:10, 10, 1)) # TRUE
#' has_dim(matrix(1:10, 1, 10)) # TRUE
#' 
#' has_dim(1) # FALSE
#' lacks_dim(1) # TRUE
#' has_dim(1:10) # FALSE
#' has_dimension("dimension") # FALSE
has_dimension <- function(x) {
  if (!is.null(dim(x))) TRUE else FALSE
}

has_dim <- function(x) {
  has_dimension(x)
}

lacks_dimension <- function(x) {
  !has_dimension(x)
}

lacks_dim <- function(x) {
  !has_dimension(x)
}
