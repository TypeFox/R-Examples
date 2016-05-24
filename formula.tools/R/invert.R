#' invert
#' 
#' Invert the operators in an object, usually a formula or expression
#' 
#' @param x function for invert 
#' @param ... additional arguments passed other functions
#' 
#' \code{invert} is a S4 generic method for inverting relational
#' operators, i.e. 
# changing \code{>} to \code{<=} and \code{%in%} to 
# \code{%!in%} etc.
#' 
#' functions prefixed with a \code{.} are not exported and should probably not 
#' be called directly
#' 
#' @return 
#'   The operand is returned with the relational operators inverted.
#' 
#' @seealso 
#'   \code{\link{op}}, \code{\link{op.type}} 
#' 
#' @author Christopher Brown
#' 
#' @examples
#'   invert( quote( A >  5 ) )
#'   invert( quote( A >= 5 ) )
#'   invert( quote( A <  5 ) )
#'   invert( quote( A <= 5 ) )
#'   invert( quote( A == 5 ) )
#'   invert( quote( A != 5 ) )
#'   invert( quote( A %in% lettters[1:5] ) )
#'   invert( quote( A %!in% letters[1:5] ) )
#' 
#' @rdname invert-methods
#' @export
#' @name invert 

# if( ! isGeneric( 'invert' ) ) {
  setGeneric( 'invert', function(x, ...) standardGeneric( 'invert' ) )
# }  


# Before declaring a new generic function we check to see if it exists.
#  package::hash already defines an invert generic, so this is not 
#  necessary.


#' @rdname invert
#' @aliases .invert.single
.invert.single <- 
  function(x) { 

    o <- as.character(op(x)) 
    
    if ( o %in% operators( types="relational" ) ) {
      op(x) <- as.name( .Options$operators[[o]][['inverse']] )
    } else  {
      warning( "No inverse found for op:", op(x) )    
    }

    return(x)
}


#' @rdname invert-methods
#' @aliases invert,call-method 
setMethod( 'invert', 'call', .invert.single ) 


#' Invert multiple elements of a multiple element object
#' 
#' @rdname invert
#' @param x object to invert from
#' @seealso .invert.single
.invert.plural <- 
  function(x) {

    for( i in 1:length(x) )
      x[[i]] <- invert( x[[i]] )
    
    return(x)

  }   


#' @rdname invert-methods
#' @aliases invert,expression-method
setMethod( 'invert', 'expression', .invert.plural )

