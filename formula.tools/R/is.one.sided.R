# ---------------------------------------------------------------------
# is.one.sided
#  This function identifies unary operators/functions.  
#  - AFAIK, the only unary operators are the one-sided version of '~' 
#    and the logical negation operator '!'
#
# ---------------------------------------------------------------------

#' Determine if an object is one- or two- sided.
#'  
#' Test whether a object (typically formula, call or expression) is one-, e.g. 
#' \code{ ~ RHS } or two-sided, e.g. \code{LHS = RHS} 
#'
#' @param x object to test for one-sidedness.  
#' @param ... arguments passed to called functions
#'
#' \code{.is.one.sided} and \code{.is.one.sided.plural} are not exported and 
#' should not be relied upon. This may change in the near future.
#' 
#' @return logical; whether \code{x} is one sided or 
#' @aliases .is.one.sided
#' @rdname is.one.sided 

.is.one.sided <- function(x, ...) 
    is.name(x[[1]])  && 
    deparse(x[[1]]) %in% c( '~', '!') && 
     length(x) == 2 


#' @docType methods 
#' @rdname is.one.sided
#' @export is.one.sided
setGeneric( 'is.one.sided', function(x, ...) standardGeneric( 'is.one.sided' ) )


#' @rdname is.one.sided
#' @aliases is.one.sided,formula-method
setMethod( 'is.one.sided', 'formula', .is.one.sided )

#' @rdname is.one.sided
#' @aliases is.one.sided,call-method
setMethod( 'is.one.sided', 'call', .is.one.sided )

#' @rdname is.one.sided
#' @aliases is.one.sided,<--method
setMethod( 'is.one.sided', '<-', .is.one.sided )


# PLURAL METHODS
#' @rdname is.one.sided
#' @aliases .is.one.sided.plural
.is.one.sided.plural <- function(x,...) sapply(x, is.one.sided) 

#' @rdname is.one.sided
#' @aliases is.one.sided,expression-method
setMethod( 'is.one.sided', 'expression', .is.one.sided.plural ) 

#' @rdname is.one.sided
#' @aliases is.one.sided,list-method
setMethod( 'is.one.sided', 'list', .is.one.sided.plural )

#' @rdname is.one.sided
#' @aliases is.one.sided,ANY-method
setMethod( 'is.one.sided', 'ANY', 
  function(x, ...) {
    warning( "'is.one.sided' is not defined for object of class: ", class(x) )
    return(NA)
  }
)
