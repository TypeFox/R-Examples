#' Define boundary for pattern matching
#' 
#' Sets boundary type matching 
#' 
#' @param object target or pattern for search 
#' @param type character; one of parial (default), full, word, sentence, line,
#'                        starts_with, ends_with
#' 
#' Sets the options for matching at specified boundaries.
#' 
#' Boundaries may also be supplied to by pattern types; regex, fixed, coll, ,,,
#' When declared explicitly, these take precedent.  
#' 
#' Except when the boundad
#' 
#' @seealso 
#'   # -tk
#' 
#' @examples
#'   # -tk
#' 
#' @rdname boundary
#' @export

boundary <- function( 
    object
  , type = c('partial','full','word','sentence','line', 'starts_with', 'ends_with') 
) {
  type = match.arg(type)
  object %<>% .set_boundary(type)  
}  

#' @rdname boundary
#' @export
full        <- function(object) boundary(object,'full')

#' @rdname boundary
#' @export
partial     <- function(object) boundary(object,'partial')

#' @rdname boundary
#' @export
word        <- function(object) boundary(object,'word')

#' @rdname boundary
#' @export
sentence    <- function(object) boundary(object,'sentence')


#' @rdname boundary
#' @export
startswith <- function(object) boundary(object,'starts_with')

#' @rdname boundary
#' @export
endsqwith   <- function(object) boundary(object,'ends_with')



# --------------------------------------------------------------
# UTILITIES
# --------------------------------------------------------------

.set_boundary <- function(object, which = c('full','partial','word','sentence','line', 'atarts_with', 'ends_with') ) {
  
  value <- match.arg(which)
  attr( object, "boundary" ) <- value 
  
  return(object) 
  
}