#' @rdname is_real
#' @export
is_imaginary <- function(x, tol = 100 * .Machine$double.eps,
  .xname = get_name_in_parent(x))
{
  if(!is.numeric(x))
  {
    x <- coerce_to(x, "complex", .xname)
  }
  call_and_name(
    function(x)
    {
      ok <- abs(Re(x)) < tol
      set_cause(ok, "real")
    },
    x
  )
}

#' Is the input real/imaginary?
#'
#' Checks to see if the input is real or imaginary.
#'
#' @param x Input to check.
#' @param tol Imaginary/real components smaller than \code{tol} are not 
#' considered.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{TRUE} if the input has imaginary component equal to zero.
#' The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[base]{complex}}
#' @examples
#' (x <- with(expand.grid(re = -1:1, im = -1:1), re + im * 1i))
#' is_real(x)
#' is_imaginary(x)
#' 
#' # By default, very small imaginary/real components are ignored.
#' x <- .Machine$double.eps * (1 + 1i)
#' is_real(x)
#' is_real(x, 0)
#' is_imaginary(x)
#' is_imaginary(x, 0)
#' # numbers with both a real and imaginary component return FALSE
#' # (since they are neither purely real nor purely imaginary)
#' cmplx <- 1 + 1i
#' is_real(cmplx)
#' is_imaginary(cmplx)
#' assert_all_are_real(1:10)
#' assert_all_are_real(1:10 + 0i)
#' assert_any_are_real(c(1i, 0))
#' assert_all_are_imaginary(1:10 * 1i)
#' assert_any_are_imaginary(c(1i, 0))
#' assertive.base::dont_stop(assert_all_are_real(x))
#' assertive.base::dont_stop(assert_all_are_imaginary(x))
#' @export
is_real <- function(x, tol = 100 * .Machine$double.eps,
  .xname = get_name_in_parent(x))
{
  if(!is.numeric(x))
  {
    x <- coerce_to(x, "complex", .xname)
  }
  call_and_name(
    function(x)
    {
      ok <- abs(Im(x)) < tol
      set_cause(ok, "imaginary")
    },
    x
  )
}
