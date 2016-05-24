#' @title Is matrix
#' @description
#' \code{is_matrix} tests if an object is a matrix \cr
#' \code{is_numeric_matrix} tests if an object is a numeric matrix \cr
#' \code{is_string_matrix} tests if an object is a string matrix \cr
#' \code{is_logical_matrix} tests if an object is a logical matrix \cr
#' \code{is_not_matrix} tests if an object is not a matrix
#' 
#' @param x an R object
#' @name is_matrix
#' @aliases is_matrix is_numeric_matrix is_string_matrix is_logical_matrix
#' is_not_matrix
#' @export is_matrix is_numeric_matrix is_string_matrix is_logical_matrix
#' is_not_matrix
#' @examples
#' A = matrix(1:10, 5, 2)
#' B = matrix(letters[1:10], 5, 2)
#' C = 1:10
#' 
#' is_matrix(A) # TRUE
#' is_matrix(C) # FALSE
#' is_not_matrix(C) # TRUE
#' 
#' is_numeric_matrix(A) # TRUE
#' is_numeric_matrix(B) # FALSE
#' 
#' is_string_matrix(A) # FALSE
#' is_string_matrix(B) # TRUE
NULL

is_matrix <- function(x) {
  is.matrix(x)
}

is_numeric_matrix <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  is.numeric(x)
}

is_string_matrix <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  is.character(x)
}

is_logical_matrix <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  is.logical(x)
}

is_not_matrix <- function(x) {
  !is_matrix(x)
}
