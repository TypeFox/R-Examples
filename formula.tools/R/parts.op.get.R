# -----------------------------------------------------------------------------
# op
#   extract and manipulate the operator of a call, expression or rule
#   generics are provided in the expressions package
# -----------------------------------------------------------------------------

#' @rdname formula.parts
#' @aliases op
#' @export op

setGeneric( 'op', function(x) standardGeneric( 'op' ) )

#' @rdname formula.parts
#' @aliases op,formula-method
setMethod( 'op', 'formula', function(x) x[[1]] )

#' @rdname formula.parts
#' @aliases op,call-method
setMethod( 'op', 'call', function(x) x[[1]] )

#' @rdname formula.parts
#' @aliases op,name-method
setMethod( 'op', 'name', function(x) if( as.character(x) %in% operators( "ALL" ) ) return(x) )

#' @rdname formula.parts
#' @aliases op,expression-method
setMethod( 'op', 'expression', 
  function(x) {
    ret <- vector( 'expression', length(x) )
    for( i in 1:length(x) ) {
      o <- op( x[[i]] ) 
      if( ! is.null(op) ) ret[[i]] <- o
    }
    ret
  }
)

#' @rdname formula.parts
#' @aliases op,list-method
setMethod( 'op', 'list', function(x) lapply(x,op) )


#' @rdname formula.parts
#' @aliases op,<--method
setMethod( 'op', '<-', function(x) x[[1]] ) 

