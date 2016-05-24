#' @rdname is_in_range
#' @export
is_in_closed_range <- function(x, lower = -Inf, upper = Inf, 
  .xname = get_name_in_parent(x))
{
  is_in_range(x, lower, upper, FALSE, FALSE, .xname = .xname)
}

#' @rdname is_in_range
#' @export
is_in_left_open_range <- function(x, lower = -Inf, upper = Inf, 
  .xname = get_name_in_parent(x))
{
  is_in_range(x, lower, upper, TRUE, FALSE, .xname = .xname)
}

#' @rdname is_in_range
#' @export
is_in_open_range <- function(x, lower = -Inf, upper = Inf, 
  .xname = get_name_in_parent(x))
{
  is_in_range(x, lower, upper, TRUE, TRUE, .xname = .xname)
}

#' Is the input in range?
#'
#' Checks to see if the input is within an numeric interval.
#'
#' @param x Input to check.
#' @param lower Lower bound for the interval.
#' @param upper Upper bound for the interval.
#' @param lower_is_strict If \code{TRUE}, the lower bound is open (strict) 
#' otherwise it is closed.
#' @param upper_is_strict If \code{TRUE}, the upper bound is open (strict)
#' otherwise it is closed.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @note \code{is_in_range} provides the most flexibility in determining
#' if values are within a numeric interval.  The other functions restrict
#' the input arguments for convience in common cases.  For example,
#' \code{is_percentage} forces the interval to be from 0 to 100.
#' The function is not vectorized by the \code{lower_is_strict} and
#' \code{upper_is_strict} for speed (these are assumed to be scalar logical
#' values).
#' @return The \code{is_*} functions return \code{TRUE} if the input is 
#' within an interval.  The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @examples
#' assert_all_are_positive(1:10)
#' assert_all_are_non_negative(0:10)
#' assert_any_are_positive(c(-1, 1))
#' assert_all_are_percentages(c(0, 50, 100))
#' assert_all_are_proportions(c(0, 0.5, 1))
#' assert_all_are_in_left_open_range(1 + .Machine$double.eps, lower = 1)
#' @importFrom assertive.base use_first
#' @export
is_in_range <- function(x, lower = -Inf, upper = Inf, lower_is_strict = FALSE, 
  upper_is_strict = FALSE, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  lower <- coerce_to(lower, "numeric")
  upper <- coerce_to(upper, "numeric")
  lower_is_strict <- coerce_to(use_first(lower_is_strict), "logical")
  upper_is_strict <- coerce_to(use_first(upper_is_strict), "logical")
  ok <- rep.int(TRUE, length(x))
  ok[is.na(x)] <- NA
  too_low <- (if(lower_is_strict) `<=` else `<`)(x, lower)
  too_high <- (if(upper_is_strict) `>=` else `>`)(x, upper)
  ok[too_low] <- FALSE                     
  ok[too_high] <- FALSE
  names(ok) <- x
  set_cause(
    ok,
    ifelse(too_low, "too low", "too high")
  )
}

#' @rdname is_in_range
#' @export
is_in_right_open_range <- function(x, lower = -Inf, upper = Inf, 
  .xname = get_name_in_parent(x))
{
  is_in_range(x, lower, upper, FALSE, TRUE, .xname = .xname)
}

#' @rdname is_in_range
#' @export
is_negative <- function(x, .xname = get_name_in_parent(x))
{
  is_in_range(x, upper = 0, upper_is_strict = TRUE, .xname = .xname)
}

#' @rdname is_in_range
#' @export
is_non_negative <- function(x, .xname = get_name_in_parent(x))
{
  is_in_range(x, 0)
}

#' @rdname is_in_range
#' @export
is_non_positive <- function(x, .xname = get_name_in_parent(x))
{
  is_in_range(x, upper = 0, .xname = .xname)
}

#' @rdname is_in_range
#' @export
is_percentage <- function(x, lower_is_strict = FALSE, upper_is_strict = FALSE, 
  .xname = get_name_in_parent(x))
{
  is_in_range(x, 0, 100, lower_is_strict, upper_is_strict, .xname = .xname)
}

#' @rdname is_in_range
#' @export
is_positive <- function(x, .xname = get_name_in_parent(x))
{
  is_in_range(x, 0, lower_is_strict = TRUE, .xname = .xname)
}

#' @rdname is_in_range
#' @export
is_proportion <- function(x, lower_is_strict = FALSE, upper_is_strict = FALSE, 
  .xname = get_name_in_parent(x))
{
  is_in_range(x, 0, 1, lower_is_strict, upper_is_strict, .xname = .xname)
}
