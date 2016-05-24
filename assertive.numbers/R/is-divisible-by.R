#' Is the input divisible by a number?
#' 
#' Checks to see if the input is divisible by some number.
#' @param x A numeric vector to divide.
#' @param n A numeric vector to divide by.
#' @param tol Differences from zero smaller than \code{tol} are not considered.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{TRUE} if the input \code{x} is divisible by \code{n}, within 
#' the specified tolerance.
#' @note \code{is_even} and \code{is_odd} are shortcuts for divisibility by two.
#' @seealso \code{is_whole_number}
#' @examples
#' is_divisible_by(1:10, 3)
#' is_divisible_by(-5:5, -2)
#' is_divisible_by(1.5:10.5, c(1.5, 3.5))
#' assert_any_are_even(1:10)
#' assertive.base::dont_stop(assert_all_are_even(1:10))
#' @export
is_divisible_by <- function(x, n, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x))
{
  if(!is.integer(x))
  {
    x <- coerce_to(x, "numeric", .xname)
  }
  call_and_name(
    function(x) 
    {
      ok <- abs(x %% n) <= tol
      set_cause(ok, "indivisible")
    }, 
    x
  ) 
}

#' @rdname is_divisible_by
#' @export
is_even <- function(x, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x))
{
  is_divisible_by(x, 2L, tol = tol, .xname = .xname)  
}

#' @rdname is_divisible_by
#' @importFrom stats setNames
#' @export
is_odd <- function(x, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x))
{
  setNames(is_divisible_by(x - 1, 2L, tol = tol, .xname = .xname), x)
}

