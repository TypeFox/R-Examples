#' Does the input have terms?
#'
#' Checks to see if the input has a terms component or attribute.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_terms} returns \code{TRUE} if the input has an 
#' element or an attribute named terms. \code{assert_has_terms} returns 
#' nothing but throws an error if \code{has_terms} is not \code{TRUE}.
#' @seealso \code{\link[stats]{terms}}.
#' @examples
#' assert_has_terms(lm(uptake ~ conc, CO2))
#' @export
has_terms <- function(x, .xname = get_name_in_parent(x))
{
  if(
    is.null(attr(x, "terms")) && 
    (is.atomic(x) || is.null(x$terms))
  )
  {
    return(
      false(
        gettext("%s has no terms component nor attribute."), 
        .xname
      )
    )
  }
  TRUE
}
