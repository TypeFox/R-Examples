#' @title Is natural
#' @description Test if is a natural number
#' @details Zero is not included in the set of natural numbers  
#' @param x an R object
#' @export is_natural is_not_natural
#' @aliases is_natural is_not_natural
#' @seealso \code{\link{is_negative}}
#' @examples
#' is_natural(1)
#' is_natural(0)
#' is_natural(seq(-2, 3))
#' is_natural(iris$Species)
#' 
#' M = matrix(seq(-3, 2), 2, 3)
#' is_natural(M)
is_natural <- function(x) {
  UseMethod("is_natural", x)
}

#' @S3method is_natural default
is_natural.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_natural factor
is_natural.factor <- function(x) {
  FALSE
}

#' @S3method is_natural numeric
is_natural.numeric <- function(x) {
  ints <- (x %% 1) == 0
  ints & (x > 0)
}

is_not_natural <- function(x) {
  !is_natural(x)
}
