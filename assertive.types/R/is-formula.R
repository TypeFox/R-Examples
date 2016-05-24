#' Is the input a formula?
#'
#' Checks to see if the input is a formula.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{is_*} functions return \code{TRUE} when the input is a 
#' formula.  The \code{assert_*} functions return nothing but throw an error 
#' if the corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link{is_environment}} and  \code{\link{is_language}}
#' @examples
#' is_one_sided_formula(~ x)
#' is_two_sided_formula(y ~ x)
#' @importFrom assertive.base is2
#' @export
is_formula <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "formula", .xname)
}

#' @rdname is_formula
#' @importFrom assertive.properties is_of_length
#' @export
is_one_sided_formula <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_formula(x, .xname)))
  {
    return(ok)
  }
  if(!(ok <- is_of_length(x, 2L, .xname)))
  {
    return(ok)
  } 
  TRUE
}

#' @rdname is_formula
#' @importFrom assertive.properties is_of_length
#' @export
is_two_sided_formula <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_formula(x, .xname)))
  {
    return(ok)
  }
  if(!(ok <- is_of_length(x, 3L, .xname)))
  {
    return(ok)
  } 
  TRUE
}
