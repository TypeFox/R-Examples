#' Reject
#'
#' \code{Reject()} is the opposite of Filter.
#' Reject applies the negation of the unary predicate function f to each
#' element of x, coercing to logical if necessary, and returns the subset
#' of x for which this gives true. Note that possible NA values are currently
#' always taken as false; control over NA handling may be added in the future.
#'
#' @param f a predicate function.
#' @param x a vector.
#' @return x filtered where f applies
#' @family predicate functionals
#' @examples
#' # Some examples
#' Filter(function(x) x < 5, 1:10)
#' Reject(function(x) x < 5, 1:10)
#' @export
Reject <- function (f, x)
{
  ind <- as.logical(unlist(lapply(x, Negate(f))))
  return(x[!is.na(ind) & ind])
}
