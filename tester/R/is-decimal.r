#' @title Is decimal
#' @description Test if is a decimal number
#' @details decimal is any number in the intervals (-1,0) and (0,1) 
#' @param x an R object
#' @export is_decimal is_not_decimal
#' @aliases is_decimal is_not_decimal
#' @seealso \code{\link{is_integer}}
#' @examples
#' is_decimal(0.01) # TRUE
#' is_decimal(-0.01) # TRUE
#' is_decimal(0) # FALSE
#' is_decimal(1) # FALSE
#' is_decimal(runif(5))
#' is_decimal(rnorm(5))
#' 
#' M = matrix(seq(-2, 2, length.out=10), 5, 2)
#' is_decimal(M)
is_decimal <- function(x) {
  UseMethod("is_decimal", x)
}

#' @S3method is_decimal default
is_decimal.default <- function(x) {
  if (mode(x) != "numeric") FALSE
}

#' @S3method is_decimal factor
is_decimal.factor <- function(x) {
  FALSE
}

#' @S3method is_decimal numeric
is_decimal.numeric <- function(x) {
  (abs(x) > 0 & abs(x) < 1) 
}

is_not_decimal <- function(x) {
  !is_decimal(x)
}


#' @title Is positive decimal
#' @description Test if is a positive decimal
#' @param x an R object
#' @export
#' @examples
#' is_positive_decimal(0.0001)
#' is_positive_decimal(-0.0003)
#' is_positive_decimal(0)
#' is_positive_decimal(pi)
#' is_positive_decimal(-exp(1))
is_positive_decimal <- function(x) {
  UseMethod("is_positive_decimal", x)
}

#' @S3method is_positive_decimal default
is_positive_decimal.default <- function(x) {
  if (is_positive(x) & is_decimal(x)) TRUE else FALSE
}


#' @title Is negative decimal
#' @description Test if is a negative decimal
#' @param x an R object
#' @export
#' @examples
#' is_negative_decimal(0.0001)
#' is_negative_decimal(-0.0003)
#' is_negative_decimal(0)
#' is_negative_decimal(pi)
#' is_negative_decimal(-exp(1))
is_negative_decimal <- function(x) {
  UseMethod("is_negative_decimal", x)
}

#' @S3method is_negative_decimal default
is_negative_decimal.default <- function(x) {
  if (is_negative(x) & is_decimal(x)) TRUE else FALSE
}
