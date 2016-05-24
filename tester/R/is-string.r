#' @title Is string
#' @description Tests if an object is a character string \cr
#' \code{is_not_string()} tests the opposite condition
#' 
#' @param x an R object
#' @aliases is_string is_not_string
#' @export is_string is_not_string
#' @examples
#' is_string("test_me") # TRUE
#' 
#' is_string(1:10) # FALSE
is_string <- function(x) {
  if (is.character(x)) TRUE else FALSE
}

is_not_string <- function(x) {
  !is_string(x)
}
