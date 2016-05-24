#' @title Same Type
#' 
#' @description
#' \code{same_type()} tests if two objects have the same type \cr
#' \code{different_type()} tests if two objects have different type
#' 
#' @param x an R object
#' @param y an R object
#' @aliases same_type different_type
#' @export same_type different_type
#' @examples
#' same_type(letters[1:3], "class") # TRUE
#' same_type(1:3, "class") # FALSE
#' 
#' different_type(1, 1L) # TRUE
#' different_type(1, 1.0) # FALSE
same_type <- function(x, y) {
  (typeof(x) == typeof(y))
}

different_type <- function(x, y) {
  !same_type(x, y)
}
