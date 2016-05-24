#' coll 
#' 
#' Creates or modifies the search type to use coll matching 
#' 
#' @param object to make specification 
#' @param ... additional arguments passes to \code{\link{pattern}}
#' 
#' \code{coll}
#' 
#' @seealso 
#'   \code{\link{pattern}}
#' 
#' @examples
#'   pat <- coll("a") 
#'   detect( c('alpha','beta'), pat )
#'   

#' @rdname coll
#' @export
   coll <- function( object, ... )  UseMethod('coll')
  

#' @rdname coll
#' @export
   coll.default  <- function( object, ... ) pattern( as.character(object), 'coll', ... )


#' @rdname coll
#' @export
   coll.character  <- function( object, ... ) pattern( object, 'coll', ... )


#' @rdname coll
#' @export
#  Note: an alternative is to clobber the options
   coll.SearchableOrPattern  <- function( object, ... ) {
     object@type = 'coll'
     object@options = list(...)
     return(object)
   }
