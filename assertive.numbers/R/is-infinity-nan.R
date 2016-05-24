#' Are the inputs (in)finite?
#'
#' Checks to see if the inputs are (in)finite.
#'
#' @param x Input to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_finite} wraps \code{is.finite}, showing the 
#' names of the inputs in the answer. \code{is_infinite} works 
#' likewise for \code{is.infinite}.  The \code{assert_*} functions 
#' return nothing but throw an error if the corresponding 
#' \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{is.finite}}
#' @examples
#' x <- c(0, Inf, -Inf, NA, NaN)
#' is_finite(x)
#' is_infinite(x)
#' is_positive_infinity(x)
#' is_negative_infinity(x)
#' assert_all_are_finite(1:10)
#' assert_any_are_finite(c(1, Inf))
#' assert_all_are_infinite(c(Inf, -Inf))
#' assertive.base::dont_stop(assert_all_are_finite(c(0, Inf, -Inf, NA, NaN)))
#' @export
is_finite <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  call_and_name(
    function(x)
    {
      ok <- is.finite(x)
      set_cause(
        ok, 
        ifelse(
          is.infinite(x),
          "infinite",
          ifelse(is.nan(x), "not a number", "missing")
        )
      )
    }, 
    x
  )
}

#' @rdname is_finite
#' @export
is_infinite <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  call_and_name(
    function(x)
    {
      ok <- is.infinite(x)
      set_cause(
        ok, 
        ifelse(
          is.finite(x),
          "finite",
          ifelse(is.nan(x), "not a number", "missing")
        )
      )
    }, 
    x
  )
}

#' Is the input (not) NaN?
#'
#' Checks to see if the input is a number that is(n't) NaN.
#'
#' @param x Input to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_nan} wraps \code{is.nan}, coercing the input to
#' numeric if necessary.  \code{is_not_nan} works similarly, but returns
#' the negation.  The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns
#' \code{FALSE}.
#' @seealso \code{\link[base]{is.nan}}
#' @examples
#' x <- c(0, NaN, NA)
#' is_nan(x)
#' is_not_nan(x)
#' assert_all_are_not_nan(1:10)
#' assert_any_are_not_nan(x)
#' assertive.base::dont_stop(assert_all_are_not_nan(x))
#' @export
is_nan <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  call_and_name(
    function(x)
    {
      ok <- is.nan(x)
      set_cause(ok, "a number")
    }, 
    x
  )
}

#' @rdname is_nan
#' @export
is_not_nan <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  call_and_name(
    function(x)
    {
      ok <- !is.nan(x)
      set_cause(ok, "not a number")
    }, 
    x
  )
}

#' @rdname is_finite
#' @export
is_negative_infinity <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  call_and_name(
    function(x)
    {
      ok <- is.infinite(x) & x < 0
      set_cause(
        ok, 
        ifelse(
          is.finite(x),
          "finite",
          ifelse(
            is.nan(x), 
            "not a number",
            ifelse(is.na(x), "missing", "positive inf")
          )
        )
      )
    }, 
    x
  )
}

#' @rdname is_finite
#' @export
is_positive_infinity <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  call_and_name(
    function(x)
    {
      ok <- is.infinite(x) & x > 0
      set_cause(
        ok, 
        ifelse(
          is.finite(x),
          "finite",
          ifelse(
            is.nan(x), 
            "not a number",
            ifelse(is.na(x), "missing", "negative inf")
          )
        )
      )
    }, 
    x
  )
}
