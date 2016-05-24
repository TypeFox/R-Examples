#' @title Same Class
#' 
#' @description
#' \code{same_class()} tests if two objects have the same class \cr
#' \code{different_class()} tests if two objects have different class
#' 
#' @param x an R object
#' @param y an R object
#' @aliases same_class different_class
#' @export same_class different_class
#' @examples
#' same_class(letters[1:3], "class") # TRUE
#' same_class(1:3, "class") # FALSE
same_class <- function(x, y) {
  identical(class(x), class(y))
}

different_class <- function(x, y) {
  !same_class(x, y)
}
