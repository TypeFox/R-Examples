# -----------------------------------------------------------------------------
# op.type
#   extract and manipulate the operator type of a call, expression, or rule
#
#   n.b.
#     - No replacement methods available for op.type.
# -----------------------------------------------------------------------------

# OP.TYPE

#' Get the operator type used in an call, formula, expression, etc.
#' 
#' @param x object from which to extract the operator type
#' 
#' @return a character vector of the operator type(s) 
#' 
#' @seealso \code{\link{op}}, \code{\link[operator.tools]{operator.type}}
#' @docType methods
#' @rdname op.type
#' @aliases op.type-methods
#' @export op.type

setGeneric( 'op.type', function(x) standardGeneric( 'op.type' ) )


# SINGULAR METHODS

#' @rdname op.type
#' @aliases op.type,call-method
setMethod( 'op.type', 'call' , function(x) operator.type( (op(x) )  ) )

#' @rdname op.type
#' @aliases op.type,formula-method
setMethod( 'op.type', 'formula' , function(x) operator.type( op(x) ) )

#' @rdname op.type
#' @aliases op.type,`<-`-method
setMethod( 'op.type', '<-', function(x) operator.type( op(x) ) ) 

#' @rdname op.type
#' @aliases op.type,ANY-method
setMethod( 'op.type', 'ANY', function(x) operator.type( op(x) ) )



# PLURAL METHODS
#' @rdname op.type
#' @aliases op.type,expression-method
setMethod(  'op.type', 'expression', 
  function(x) lapply( x, function(x) operator.type( op(x) )  )  
)

#' @rdname op.type
#' @aliases op.type,list-method
setMethod( 'op.type', 'list',
  function(x) lapply( x, function(x) operator.type( op(x) )  )  
)
                           