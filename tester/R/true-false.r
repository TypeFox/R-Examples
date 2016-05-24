#' @title If TRUE or FALSE
#' @description
#' \code{is_TRUE} and \code{is_true} tests if x is TRUE \cr
#' \code{is_FALSE} and \code{is_false} tests if x is FALSE \cr
#' \code{true_or_false} returns whether the condition is true or false
#' 
#' @param x an R object
#' @name is_TRUE
#' @aliases is_TRUE is_FALSE is_true is_false true_or_false
#' @export is_TRUE is_FALSE is_true is_false true_or_false
#' @examples
#' is_true(TRUE)
#' is_true(FALSE)
#' is_false(TRUE)
#' is_false(FALSE)
#' true_or_false(TRUE)
#' true_or_false(FALSE)
#' 
#' is_true(1) # FLASE
#' is_false("FALSE") # FALSE
NULL

is_TRUE <- function(x) {
  if (is.logical(x)) {
    if (is.na(x)) {
      FALSE
    } else {
      if (x == TRUE) TRUE else FALSE      
    }
  } else FALSE
}

is_FALSE <- function(x) {
  if (is.logical(x)) {
    if (is.na(x)) {
      FALSE
    } else {
      if (x == FALSE) TRUE else FALSE
    }
  } else FALSE
}

is_true <- function(x) is_TRUE(x)

is_false <- function(x) is_FALSE(x)

true_or_false <- function(x) {
  if (is.logical(x)) x else FALSE
}
