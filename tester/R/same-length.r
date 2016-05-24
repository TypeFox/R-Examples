#' @title Same Length
#' 
#' @description 
#' \code{same_length()} tests if two objects have same length \cr
#' \code{different_length()} tests if two objects have different length
#' 
#' @param x a matrix
#' @param y a matrix
#' @aliases same_length different_length
#' @export same_length different_length
#' @examples
#' same_length(1:10, letters[11:20]) # TRUE
#' same_length(1:10, letters[11:19]) # FALSE
#' 
#' a = matrix(1:15, 5, 3)
#' same_length(a, a) # TRUE
#' same_length(a, t(a)) # TRUE
#' 
#' different_length(t(a), a) # FALSE
#' different_length(1:10, a) # TRUE
#' different_length(a, "a") # TRUE
same_length <- function(x, y) {
  (length(x) == length(y))
}

different_length <- function(x, y) {
  !same_length(x, y)
}
