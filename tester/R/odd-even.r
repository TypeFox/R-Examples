#' @title Is even
#' @description Test if an object is an even number \cr
#' \code{is_not_even} tests the opposite condition
#' 
#' @param x an R object
#' @export is_even is_not_even
#' @aliases is_even is_not_even
#' @seealso \code{\link{is_odd}}
#' @examples
#' is_even(2)
#' is_even(1)
#' is_even(seq(-5, 5))
#' 
#' is_even(iris$Species)
#' is_even(iris)
#' is_even(list(1, 0, -1, iris))
#' 
#' set.seed(999)
#' M = matrix(1:12, 4, 3)
#' is_even(M)
is_even <- function(x) {
  UseMethod("is_even", x)
}

#' @S3method is_even default
is_even.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_even numeric
is_even.numeric <- function(x) {
  (x %% 2) == 0
}

#' @S3method is_even matrix
is_even.matrix <- function(x) {
  (x %% 2) == 0
}

#' @S3method is_even factor
is_even.factor <- function(x) {
  FALSE
}

is_not_even <- function(x) {
  !is_even(x)
}

#' @title Is even
#' @description Test if an object is an even number \cr
#' \code{is_not_odd} tests the opposite condition
#' @param x an R object
#' @export is_odd is_not_odd
#' @aliases is_odd is_not_odd
#' @seealso \code{\link{is_even}}
#' @examples
#' is_odd(2)
#' is_odd(1)
#' is_odd(seq(-5, 5))
#' 
#' is_odd(iris$Species)
#' is_odd(iris)
#' is_odd(list(1, 0, -1, iris))
#' 
#' set.seed(999)
#' M = matrix(1:12, 4, 3)
#' is_odd(M)
is_odd <- function(x) {
  UseMethod("is_odd", x)
}

#' @S3method is_odd default
is_odd.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_odd numeric
is_odd.numeric <- function(x) {
  (x %% 2) != 0
}

#' @S3method is_odd matrix
is_odd.matrix <- function(x) {
  (x %% 2) != 0
}

#' @S3method is_odd factor
is_odd.factor <- function(x) {
  FALSE
}

is_not_odd <- function(x) {
  !is_odd(x)
}
