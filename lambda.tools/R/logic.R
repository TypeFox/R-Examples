# :vim set filetype=R

#' Check if an object is a scalar
#'
#' This function checks if an object is a scalar.
#'
#' @section Usage:
#' is.scalar(x)
#'
#' @section Details:
#' This function checks to determine if an object \code{x} is a scalar,
#' i.e. the length of the object is equal to one.
#'
#' @name is.scalar
#' @param x An object
#' @return A logical value that indicates if the input is of length one
#'
#' @examples
#' is.scalar(10)
#' 
#' is.scalar(1:10)
is.scalar(x) %when% { length(x) == 1 } %as% TRUE
is.scalar(x) %as% FALSE

#' Conditionally apply a function to an argument
#'
#' This function conditionally applies a function to an argument given
#' a logical condition.
#'
#' @section Usage:
#' onlyif(condition, fn, x)
#'
#' @section Details:
#' This function can be used to apply a function to a vector containing
#' elements that lie outside the valid domain of \code{fn}. 
#' The function \code{onlyif} differs from \code{ifelse} in the sense that
#' it is not vectorized and a closure can be used. For example,
#'
#' \code{ifelse(length(x) < 10, pad(x, 10 - length(x)), x)}
#'
#' yields the wrong result due to the length of the first argument. The
#' \code{onlyif} function is designed for these situations.
#'
#' \code{onlyif(length(x) < 10, function(x) pad(x, 10 - length(x)), x)}.
#'
#' Note that a closure is only required if an expression cannot be evaluated
#' under both a TRUE or FALSE scenario.
#'
#' The alternative would be to use a conditional block, which can result
#' in improperly scoped code if one is careless.
#'
#' @note
#' The interface for this function is experimental. I'm looking for a way
#' to preserve unevaluated expressions. Until then, I don't recommend
#' using the function.
#'
#' @name onlyif
#' @param condition Logical statement used to conditionally apply fn to x
#' @param fn A function to apply to x
#' @param expr An expression 
#' @param x An object
#' @return Either \code{expr} if \code{condition} is true, otherwise \code{x}.
#' @seealso \code{\link{use_default}}
#'
#' @examples
#' x <- 1:5
#' onlyif(length(x) < 10, pad(x, 10 - length(x)), x)
#' onlyif(length(x) < 10, function(x) pad(x, 10 - length(x)), x)
#'
#' # This returns x
#' x <- 1:20
#' onlyif(length(x) < 10, function(x) pad(x, 10 - length(x)), x)
onlyif(condition, fn, x) %::% logical : Function : . : .
onlyif(TRUE, fn, x) %as% fn(x)
onlyif(FALSE, fn, x) %as% x

onlyif(condition, expr, x) %::% logical : . : . : .
onlyif(TRUE, expr, x) %as% expr
onlyif(FALSE, expr, x) %as% x

#' Apply a default value whenever a variable is not well-formed
#'
#' This function provides a functional approach for a specific use case
#' of conditional expressions: that of applying default values when
#' a variable is not well-formed. In this context, well-formedness is
#' considered to be any scalar value that is not NA.
#' By encapsulating this behavior in a
#' function, referential transparency is preserved.
#'
#' @name use_default
#' @param x a scalar variable
#' @param default the value to replace empty, NULL, or NA
#'
#' @return A well-formed value, either the original value or the default
#' if x is not well-formed.
#' @seealso \code{\link{onlyif}}
#'
#' @examples
#' x <- c(1, 2, 3, NA, NA)
#' map(x, function(y) use_default(y, 0))
use_default(EMPTY, default) %as% default
use_default(NULL, default) %as% default
use_default(NA, default) %as% default
use_default(x, default) %when% { is.scalar(x) } %as% x
