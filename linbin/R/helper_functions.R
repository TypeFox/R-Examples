#' Full Name Matching
#' 
#' Recursive wrapper for \code{grep} returning only full matches to elements of a character vector.
#' 
#' @param patterns Character vector containing regular expressions to be matched.
#' @param x Character vector where matches are sought, or an object which can be coerced by \code{as.character} to a character vector.
#' @param ... Additional arguments to \code{\link{grep}}.
#' @param simplify If \code{TRUE}, the result is simplified from a list to a vector or matrix if appropriate.
#' @return List, matrix, or vector of the indices of the elements of \code{x} that yielded a match to each element of \code{patterns}.
#' @keywords internal
rgrep_exact <- function(patterns, x, ..., simplify = FALSE) {
  sapply(patterns, function(pattern) grep(paste0("^", pattern, "$"), as.character(x), ...), simplify = simplify, USE.NAMES = FALSE)
}

#' Not Integers
#' 
#' Tests whether a vector is empty, non-numeric, or contains any non-integer numbers.
#' 
#' @param x An R object.
#' @return \code{FALSE} if \code{x} is a numeric vector containing only whole numbers, \code{TRUE} otherwise.
#' @keywords internal
is_not_integer <- function(x) {
  return(!is.numeric(x) || length(x) == 0 || any((x %% 1) != 0))
}

#' Generalized If Else
#' 
#' Returns different values depending on whether a test is \code{TRUE} or \code{FALSE}.
#' 
#' @param test An object which can be coerced to a logical value.
#' @param yes Value returned if \code{test} is \code{TRUE}.
#' @param no Value returned if \code{test} is \code{FALSE}.
#' @seealso \code{\link{ifelse}}
#' @keywords internal
if_else <- function(test, yes, no) {
  if (test) {
    return(yes)
  } else {
    return(no)
  }
}
