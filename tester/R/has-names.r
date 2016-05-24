#' @title Has or lacks names?
#' @description 
#' \code{has_names} tests if an object has names \cr
#' \code{lacks_names} tests if an object lacks names
#' 
#' @param x an R object
#' @aliases has_names lacks_names
#' @export has_names lacks_names
#' @seealso \code{\link{has_rownames}}
#' @examples
#' set.seed(1)
#' x <- y <- runif(10)
#' names(x) = letters[1:10]
#' 
#' has_names(x) # TRUE
#' has_names(y) # FALSE
#' 
#' lacks_names(x) # FALSE
#' lacks_names(y) # TRUE
has_names <- function(x) {
  if (!is.null(names(x))) TRUE else FALSE
}

lacks_names <- function(x) {
  !has_names(x)
}


#' @title Has or lacks row/column names?
#' @description 
#' \code{has_rownames} tests if an object has row names \cr
#' \code{has_colnames} tests if an object has column names \cr
#' \code{has_dimnames} tests if an object has dimnames \cr
#' \code{lacks_rownames} tests if an object lacks row names \cr
#' \code{lacks_colnames} tests if an object lacks column names \cr
#' \code{lacks_dimnames} tests if an object lacks dimnames \cr
#' 
#' @param x an R object
#' @aliases has_rownames has_colnames has_dimnames 
#' lacks_rownames lacks_colnames lacks_dimnames
#' @export has_rownames has_colnames has_dimnames 
#' lacks_rownames lacks_colnames lacks_dimnames
#' @seealso \code{\link{has_names}}
#' @examples
#' has_rownames(iris) # TRUE
#' has_colnames(iris) # TRUE
#' 
#' lacks_rownames(letters[1:10]) # TRUE
#' lacks_colnames(letters[1:10]) # TRUE
#' 
#' A = matrix(1:10)
#' has_dimnames(A) # FALSE
#' lacks_dimnames(A) # TRUE
has_rownames <- function(x) {
  if (!is.null(rownames(x))) TRUE else FALSE
}

has_colnames <- function(x) {
  if (!is.null(colnames(x))) TRUE else FALSE
}

has_dimnames <- function(x) {
  if (!is.null(dimnames(x))) TRUE else FALSE
}

lacks_rownames <- function(x) {
  !has_rownames(x)
}

lacks_colnames <- function(x) {
  !has_colnames(x)
}

lacks_dimnames <- function(x) {
  !has_dimnames(x)
}
