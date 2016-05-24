#' Is the input a whole number?
#'
#' Checks that the (probably floating point) input is a whole number.
#' 
#' @param x Input to check.
#' @param tol Differences smaller than \code{tol} are not considered.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @note The term whole number is used to distinguish from integer in
#' that the input \code{x} need not have type \code{integer}.  In fact
#' it is expected that \code{x} will be \code{numeric}.
#' @return \code{TRUE} if the input is a whole number.
#' @seealso \code{is_divisible_by}
#' @examples
#' # 1, plus or minus a very small number
#' x <- 1 + c(0, .Machine$double.eps, -.Machine$double.neg.eps)
#' # By default, you get a bit of tolerance for rounding errors
#' is_whole_number(x)
#' # Set the tolerance to zero for exact matching.
#' is_whole_number(x, tol = 0)
#' @export
is_whole_number <- function(x, tol = 100 * .Machine$double.eps, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "numeric", .xname)
  call_and_name(
    function(x) 
    {
      ok <- abs(x - round(x)) <= tol & !is.infinite(x)
      set_cause(ok, ifelse(is.infinite(x), "infinite", "fractional"))
    }, 
    x
  )
}
