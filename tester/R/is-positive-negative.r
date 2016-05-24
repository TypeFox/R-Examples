#' @title Is positive
#' @description Test if an object is positive
#' @param x an R object
#' @export is_positive is_not_positive
#' @aliases is_positive is_not_positive
#' @seealso \code{\link{is_negative}}
#' @examples
#' is_positive(1)
#' is_positive(0)
#' is_positive(-1)
#' is_positive(iris$Species)
#' is_positive(iris)
#' is_positive(list(1, 0, -1, iris))
#' 
#' set.seed(999)
#' M = matrix(rnorm(12), 4, 3)
#' is_positive(M)
is_positive <- function(x) {
  UseMethod("is_positive", x)
}

#' @S3method is_positive default
is_positive.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_positive numeric
is_positive.numeric <- function(x) {
  x > 0
}

#' @S3method is_positive matrix
is_positive.matrix <- function(x) {
  x > 0
}

#' @S3method is_positive factor
is_positive.factor <- function(x) {
  FALSE
}

is_not_positive <- function(x) {
  !is_positive(x)
}



#' @title Is negative
#' @description Test if an object is negative
#' @param x an R object
#' @export is_negative is_not_negative
#' @aliases is_negative is_not_negative
#' @seealso \code{\link{is_positive}}
#' @examples
#' is_negative(1)
#' is_negative(0)
#' is_negative(-1)
#' is_negative(iris$Species)
#' is_negative(iris)
#' is_negative(list(1, 0, -1, iris))
#' 
#' set.seed(999)
#' M = matrix(rnorm(12), 4, 3)
#' is_negative(M)
is_negative <- function(x) {
  UseMethod("is_negative", x)
}

#' @S3method is_negative default
is_negative.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_negative numeric
is_negative.numeric <- function(x) {
  x < 0
}

#' @S3method is_negative matrix
is_negative.matrix <- function(x) {
  x < 0
}

#' @S3method is_negative factor
is_negative.factor <- function(x) {
  FALSE
}

is_not_negative <- function(x) {
  !is_negative(x)
}

