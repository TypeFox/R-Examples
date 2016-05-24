#' regex 
#' 
#' Creates or modifies the search type to use regular expression matching 
#' 
#' @param object to make specification 
#' @param ... additional arguments passes to \code{\link{pattern}}
#' 
#' \code{regex}
#' 
#' @seealso 
#'   \code{\link{pattern}}
#' 
#' @examples
#'   pat <- regex("a.+") 
#'   detect( c('alpha','beta'), pat )
#'   

#' @rdname regex
#' @export
   regex <- function( object, ... )  UseMethod('regex')
  

#' @rdname regex
#' @export
   regex.default  <- function( object, ... ) pattern( as.character(object), 'regex', ... )


#' @rdname regex
#' @export
   regex.character  <- function( object, ... ) pattern( object, 'regex', ... )


#' @rdname regex
#' @export
#  Note: an alternative is to clobber the options
   regex.SearchableOrPattern  <- function( object, ... ) {
     object@type = 'regex'
     object@options = list(...)
     return(object)
   }
