#' @title Is vector
#' @description
#' \code{is_vector} tests if an object is a vector \cr
#' \code{is_numeric_vector} tests if an object is a numeric vector \cr
#' \code{is_string_vector} tests if an object is a string vector \cr
#' \code{is_logical_vector} tests if an object is a logical vector \cr
#' \code{is_not_vector} tests if an object is not a vector \cr
#' 
#' @param x an R object
#' @name is_vector
#' @aliases is_vector is_numeric_vector is_string_vector 
#' is_logical_vector is_not_vector
#' @export is_vector is_numeric_vector is_string_vector 
#' is_logical_vector is_not_vector
#' @examples
#' a = 1:10
#' b = letters[1:10]
#' d = matrix(1:10, 5, 2)
#' 
#' is_vector(a) # TRUE
#' is_vector(b) # TRUE
#' is_vector(d) # FALSE
#' is_not_vector(d) # TRUE
#' 
#' is_numeric_vector(a) # TRUE
#' is_numeric_vector(b) # FALSE
#' 
#' is_string_vector(a) # FALSE
#' is_string_vector(b) # TRUE
NULL

is_vector <- function(x) {
  if (is.vector(x)) TRUE else FALSE
}

is_numeric_vector <- function(x) {
  if (is.vector(x) & is.numeric(x)) TRUE else FALSE
}

is_string_vector <- function(x) {
  if (is.vector(x) & is.character(x)) TRUE else FALSE
}

is_logical_vector <- function(x) {
  if (is.vector(x) & is.logical(x)) TRUE else FALSE
}

is_not_vector <- function(x) {
  !is_vector(x)
}