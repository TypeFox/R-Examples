#' Is the input atomic/recursive/vector?
#'
#' Checks to see if the input is a type that is atomic/recursive/vector.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_atomic}, \code{is_recursive} and \code{is_vector} wrap 
#' \code{is.atomic}, \code{is.recursive} and \code{is.vector} respectively,
#' providing more information on failure.  
#' \code{is_nested} checks for recursive objects where at least one element is
#' also recursive. \code{is_non_nested} returns \code{TRUE} for atomic objects
#' and recursive objects where no elements are recursive.
#' The \code{assert_*} functions return nothing but throw an error if the 
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{is.atomic}} and \code{\link[base]{is.vector}}.
#' @examples
#' atomic_types <- list(
#'   logical(),
#'   integer(),
#'   numeric(), 
#'   complex(),
#'   character(), 
#'   raw(),
#'   matrix(), 
#'   array(),
#'   factor(),
#'   NULL
#' )
#' for(var in atomic_types) assert_is_atomic(var)
#' 
#' recursive_types <- list(
#'   list(), 
#'   expression(),
#'   data.frame(), 
#'   y ~ x,
#'   function(){},
#'   call("sin", "pi")
#' )
#' for(var in recursive_types) assert_is_recursive(var)
#' 
#' # Names are neither atomic nor recursive
#' a_name <- as.name("x")
#' is_atomic(a_name)
#' is_recursive(a_name)
#' 
#' vector_types <- c(
#'   atomic_types[1:6], 
#'   recursive_types[1:2]
#' )
#' for(var in vector_types) assert_is_vector(var)
#' 
#' # Nested objects are recursive and have at least one recursive element
#' nested_list <- list(a = 1, b = list(2:3))
#' assert_is_nested(nested_list)
#' for(elt in nested_list) assert_is_non_nested(elt)
#' @export
is_atomic <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.atomic(x))
  {
    return(false(gettext("%s is not atomic."), .xname))
  }
  TRUE
}
#' @rdname is_atomic
#' @export
is_nested <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_recursive(x, .xname)))
  {
    return(ok)
  }
  are_rec <- bapply(x, is.recursive)
  if(!any(are_rec))
  {
    return(false(gettext("%s has no recursive elements."), .xname))
  }
  TRUE
}

#' @rdname is_atomic
#' @export
is_non_nested <- function(x, .xname = get_name_in_parent(x))
{
  are_rec <- bapply(x, is.recursive)
  if(any(are_rec))
  {
    msg <- ngettext(
      sum(are_rec),
      "Element %s of %s is recursive.",
      "Elements %s of %s are recursive."
    )
    return(false(msg, toString(which(are_rec), width = 50), .xname))
  }
  TRUE
}

#' @rdname is_atomic
#' @export
is_recursive <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.recursive(x))
  {
    return(false("%s is not recursive.", .xname))
  }
  TRUE
}

#' @rdname is_atomic
#' @export
is_vector <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.vector(x)) 
  {
    return(false("%s is not a vector.", .xname))
  }
  TRUE
}                
