#' fixed 
#' 
#' Creates or modifies the search type to use fixed matching 
#' 
#' @param object to make specification 
#' @param ... additional arguments passes to \code{\link{pattern}}
#' 
#' \code{fixed}
#' 
#' @seealso 
#'   \code{\link{pattern}}
#' 
#' @examples
#'   pat <- fixed("a") 
#'   detect( c('alpha','beta'), pat )
#'   

#' @rdname fixed
#' @export
   fixed <- function( object, ... )  UseMethod('fixed')
  

#' @rdname fixed
#' @export
   fixed.default  <- function( object, ... ) pattern( as.character(object), 'fixed', ... )


#' @rdname fixed
#' @export
   fixed.character  <- function( object, ... ) pattern( object, 'fixed', ... )


#' @rdname fixed
#' @export
#  Note: an alternative is to clobber the options
   fixed.SearchableOrPattern  <- function( object, ... ) {
     object@type = 'fixed'
     object@options = list(...)
     return(object)
   }
