#' @rdname is_logical
#' @importFrom assertive.properties is_scalar
#' @export
is_a_bool <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_logical(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}

#' @rdname is_complex
#' @importFrom assertive.properties is_scalar
#' @export
is_a_complex <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_complex(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}

#' @rdname is_numeric
#' @importFrom assertive.properties is_scalar
#' @export
is_a_double <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_double(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

#' @rdname is_numeric
#' @importFrom assertive.properties is_scalar
#' @export
is_a_number <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_numeric(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

#' @rdname is_raw
#' @importFrom assertive.properties is_scalar
#' @export
is_a_raw <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_raw(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

#' @rdname is_character
#' @importFrom assertive.properties is_scalar
#' @export
is_a_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_character(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}

#' @rdname is_integer
#' @importFrom assertive.properties is_scalar
#' @export
is_an_integer <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_integer(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

#' Does the object inherit from some class?
#' 
#' Checks to see if an object is inherited from any of the specified classes.
#' @param x Any R variable.
#' @param classes A character vector of classes.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if \code{x} inherits from at least one of the classes,
#' as determined by \code{\link[base]{inherits}}.
#' @seealso \code{\link[base]{inherits}}, \code{\link[methods]{is}}, 
#' \code{\link[assertive.base]{is2}}
#' @examples
#' x <- structure(1:5, class = c("foo", "bar"))
#' assert_is_inherited_from(x, c("foo", "baz"))
#' assertive.base::dont_stop(assert_is_inherited_from(x, c("Foo", "baz")))
#' @importFrom assertive.base bapply
#' @export
is_inherited_from <- function(x, classes, .xname = get_name_in_parent(x))
{
  ok <- bapply(classes, function(class) inherits(x, class))
  if(!any(ok)) 
  {
    msg <- ngettext(
      length(classes),
      "%s does not inherit from the class %s. It has class %s.",
      "%s does not inherit from any of the classes %s. It has class %s."
    )
    return(
      false(msg, .xname, toString(classes), toString(class(x)))
    )
  }
  TRUE
}
