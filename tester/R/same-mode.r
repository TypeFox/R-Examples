#' @title Same Mode
#' 
#' @description
#' \code{same_mode()} tests if two objects have the same mode \cr
#' \code{different_mode()} tests if two objects have different mode
#' 
#' @param x an R object
#' @param y an R object
#' @aliases same_mode different_mode
#' @export same_mode different_mode
#' @examples
#' same_mode(letters[1:3], "class") # TRUE
#' same_mode(1:3, "class") # FALSE
same_mode <- function(x, y) {
  (mode(x) == mode(y))
}

different_mode <- function(x, y) {
  !same_mode(x, y)
}
