#' @title Has missing values, NA, NaN, Inf
#' @description
#' \code{has_missing} and \code{has_NA} tests if there are 
#' missing values (\code{NA}) \cr
#' \code{has_infinite} and \code{has_Inf} tests if there are 
#' infinite values (\code{Inf, -Inf}) \cr
#' \code{has_not_a_number} and \code{has_NaN} tests if there are 
#' 'Not a Number' (\code{NaN}) \cr
#' \code{has_nas} tests if there are any of the previous ones
#' 
#' @param x an R object
#' @aliases has_missing has_infinite has_not_a_number 
#' has_NA has_Inf has_NaN has_nas
#' @export has_missing has_infinite has_not_a_number 
#' has_NA has_Inf has_NaN has_nas
#' @examples
#' has_missing(1:5) # FALSE
#' has_missing(c(1, 2, 3, 4, NA)) # TRUE
#' 
#' has_infinite(c(1, 2, Inf, 1/0))
#' has_infinite(c(-Inf, "infinite"))
#' 
#' has_not_a_number(c(1, 2, 3)) # FALSE
#' has_not_a_number(c(1, 0/0, 3)) # TRUE
#' has_not_a_number(c(NaN, pi, log(1))) # TRUE
has_missing <- function(x) {
  if (sum(is.na(x)) > 0) TRUE else FALSE
}

has_infinite <- function(x) {
  if (sum(is.infinite(x)) > 0) TRUE else FALSE
}

has_not_a_number <- function(x) {
  if (sum(is.nan(x)) > 0) TRUE else FALSE
}

has_NA <- function(x) {
  has_missing(x)
}

has_Inf <- function(x) {
  has_infinite(x)
}

has_NaN <- function(x) {
  has_not_a_number(x)
}

has_nas <- function(x) {
  has_NA(x) | has_Inf(x) | has_NaN(x)
}
