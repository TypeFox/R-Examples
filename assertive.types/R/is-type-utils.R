#' Is the input relistable?
#'
#' Checks to see if the input is relistable.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_relistable} wraps \code{is.relistable}, providing more 
#' information on failure. The \code{assert_*} functions return nothing but
#' throws an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[utils]{is.relistable}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_relistable(as.relistable(list(1,2,3)))
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_relistable(list(1,2,3)))
#' @importFrom assertive.base is2
#' @export
is_relistable <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "relistable", .xname)
}

