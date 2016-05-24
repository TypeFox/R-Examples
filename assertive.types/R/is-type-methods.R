#' Is the input the name of a (formally defined) class?
#'
#' Checks to see if the input is the name of a (formally defined) class.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_class} is a vectorised wrapper for \code{isClass}.  
#' \code{assert_is_class} returns nothing but throws an error if 
#' \code{is_class} returns \code{FALSE}.
#' @seealso \code{\link[methods]{isClass}}.
#' @examples
#' assert_all_are_classes(c("lm", "numeric"))
#' @importFrom methods isClass
#' @importFrom assertive.properties is_empty
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base bapply
#' @export
is_class <- function(x, .xname = get_name_in_parent(x))
{
  if(is_empty(x)) return(logical())
  x <- coerce_to(x, "character")
  bapply(x, methods::isClass)
}
