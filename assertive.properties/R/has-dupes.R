#' Does the input have duplicates?
#'
#' Checks to see if the input has duplicates.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_duplicates} returns \code{TRUE} if\code{anyDuplicated} 
#' is \code{TRUE}.  \code{assert_has_duplicates} returns nothing but 
#' throws an error if \code{has_duplicates} is not \code{TRUE}. 
#' \code{has_no_duplicates} is the negation of \code{has_duplicates}.
#' @seealso \code{\link{anyDuplicated}}.
#' @examples 
#' x <- sample(10, 100, replace = TRUE)
#' assert_has_duplicates(x)
#' has_no_duplicates(x)
#' @export
has_duplicates <- function(x, .xname = get_name_in_parent(x))
{
  if(!anyDuplicated(x)) 
  {
    return(false(gettext("%s has no duplicates."), .xname))
  }
  TRUE
}

#' @rdname has_duplicates
#' @export
has_no_duplicates <- function(x, .xname = get_name_in_parent(x))
{
  if(anyDuplicated(x)) 
  {
    dupe_indicies <- which(duplicated(x))
    return(
      false(
        ngettext(
          length(dupe_indicies),
          "%s has a duplicate at position %s.",
          "%s has duplicates at positions %s."
        ), 
        .xname, 
        toString(dupe_indicies, width = 100)
      )
    )
  }
  TRUE
}

