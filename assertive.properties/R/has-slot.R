#' Does the S4 input have a slot?
#' 
#' Checks to see if the object is an S4 object with a particular slot.
#' @param x Input to check.  Intended to be an S4 object.
#' @param slotname A string naming a slot to check for.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_names} returns \code{TRUE} if \code{names} is 
#' non-null. 
#' @seealso \code{\link[methods]{slot}}
#' @examples 
#' setClass("numbers", representation(foo = "numeric"))
#' x <- new("numbers", foo = 1:10)
#' has_slot(x, "foo")
#' has_slot(x, "bar")
#' has_slot(1:10, "foo")
#' @importFrom methods .hasSlot
#' @importFrom methods slotNames
#' @export
has_slot <- function(x, slotname, .xname = get_name_in_parent(x))
{
  slotname <- coerce_to(use_first(slotname), "character")
  if(!isS4(x))
  {
    return(false("%s has no slots because it is not an S4 object.", .xname))
  }
  if(!.hasSlot(x, slotname))
  {
    return(
      false(
        "%s does not have a slot named %s.  The available slots are %s.", 
        .xname,
        slotname,
        toString(slotNames(x))
      )
    )
  }
  TRUE
}
