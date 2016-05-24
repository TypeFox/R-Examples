#' @title Is tabular
#' @description
#' \code{is_tabular} tests if an object has a tabular format
#' (i.e. a matrix or data frame) \cr
#' \code{is_not_tabular} tests if an object doesn't have a tabular format
#' (i.e. not matrix nor data frame) \cr
#' \code{is_numeric_tabular} tests if an object is a numeric table 
#' (i.e. a numeric matrix or data frame) \cr
#' \code{is_string_tabular} tests if an object is a string table
#' 
#' @param x an R object
#' @name is_tabular
#' @aliases is_tabular is_numeric_tabular is_string_tabular is_not_tabular
#' @export is_tabular is_numeric_tabular is_string_tabular is_not_tabular
#' @examples
#' A = matrix(1:10, 5, 2)
#' B = matrix(letters[1:10], 5, 2)
#' C = 1:10
#' 
#' is_tabular(A) # TRUE
#' is_tabular(iris) # TRUE
#' 
#' is_numeric_tabular(A) # TRUE
#' is_numeric_tabular(iris) # FALSE
#' is_numeric_dataframe(iris[,1:4]) # TRUE
NULL

is_tabular <- function(x) {
  if (is.matrix(x) | is.data.frame(x)) {
    TRUE
  } else FALSE  
}

is_numeric_tabular <- function(x) {
  if (is_numeric_matrix(x) | is_numeric_dataframe(x)) {
    TRUE
  } else FALSE  
}

is_string_tabular <- function(x) {
  if (is_string_matrix(x) | is_string_dataframe(x)) {
    TRUE
  } else FALSE  
}

is_not_tabular <- function(x) {
  !is_tabular(x)  
}
