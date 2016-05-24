#' @title Is integer
#' @description Test if a number is an integer \cr
#' Use \code{is_not_integer} to test the opposite condition
#' 
#' @param x an R object
#' @export is_integer is_not_integer
#' @aliases is_integer is_not_integer
#' @seealso \code{\link{is_natural}}
#' @examples
#' is_integer(1) # TRUE
#' is_integer(-3) # TRUE
#' is_integer(pi) # FALSE
#' is_integer(iris$Species)
#' 
#' M = matrix(seq(-3, 2), 2, 3)
#' is_integer(M)
is_integer <- function(x) {
  UseMethod("is_integer", x)
}

#' @S3method is_integer default
is_integer.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_integer factor
is_integer.factor <- function(x) {
  FALSE
}

#' @S3method is_integer numeric
is_integer.numeric <- function(x) {
  (x %% 1) == 0
}

is_not_integer <- function(x) {
  !is_integer(x)
}


#' @title Is positive integer
#' @description Test if is a positive integer
#' @param x an R object
#' @export
#' @examples
#' is_positive_integer(1) # TRUE
#' is_positive_integer(0) # FALSE
#' is_positive_integer(pi) # FALSE
#' is_positive_integer(2.2) # FALSE
#' is_positive_integer(-1) # FALSE
is_positive_integer <- function(x) {
  (is_positive(x) & is_integer(x))
}


#' @title Is negative integer
#' @description Test if is a positive integer
#' @param x an R object
#' @export
#' @examples
#' is_negative_integer(-1) # TRUE
#' is_negative_integer(1) # FALSE
#' is_negative_integer(0) # FALSE
#' is_negative_integer(pi) # FALSE
#' is_negative_integer(2.2) # FALSE
is_negative_integer <- function(x) {
  (is_negative(x) & is_integer(x))
}
