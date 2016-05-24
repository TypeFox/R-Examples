# -----------------------------------------------------------------------------
# REPLACEMENT : OP<-
# -----------------------------------------------------------------------------

#' @rdname formula.parts
#' @aliases op<-
#' @export op<-
setGeneric( 'op<-', function(x,value) standardGeneric('op<-') )

# -------------------------------------
# SINGLE: call, formula
#  - Note: if value == '~' should we eval x to return a formula?
# -------------------------------------   

#' @rdname formula.parts
#' @aliases op<-,call-method
#' @name op<-
setReplaceMethod( 'op', 'call', 
                  function( x, value ) {
                    x[[1]] <- as.name(value)
                    x
                  }
)


# EXPERIMENTAL!!!
#   Unsure of the proper behavior. Should changing of the operator for
#   a formula produce an error or should it ERROR.
# ----------------------------------------------------------------------
# METHOD: op<-,formula
#   This is a bit strange since the formula is dependent upon
#   the operator type. So if the operator is changed, we 
#   no longer have a formula, but a call object.  
#   That is, a formula appears to inherit a call.  
# ----------------------------------------------------------------------

#' @rdname formula.parts
#' @aliases op<-,formula-method
#' @name op<-
setReplaceMethod( 'op', 'formula', 
                  function( x, value ) {
                    new.op <- as.name(value) 
                    
                    # THIS CATCHES THAT WE DON"T CHANGE THE TILDE~:
                    if ( new.op == op(x) ) return(x)  
                    
                    # When we change from a tilde the operator type gets degraded.
                    if( as.character(value) %in% operators( "ALL" ) ) {
                      c <- quote( x == y )  # generic call object
                      lhs(c) <- lhs(x) 
                      op(c)  <- new.op
                      rhs(c) <- rhs(x) 
                    } else {
                      stop( value, " was not found as an operator." )
                    }
                    
                    return(c) 
                  }
)



# -------------------------------------
# LIST AND VECTORS: expression, list
# -------------------------------------

#' @rdname formula.parts
#' @aliases .replace.op.plural 
.replace.op.plural <- function( x, value ) {
  
  if( length(value) == 1  ) { 
    for( i in 1:length(x) ) op( x[[i]] ) <- as.name(value) 
    
  } else if( length(x) == length(value) ) {
    for( i in 1:length(x) ) op( x[[i]] ) <- as.name( value[[i]] )
    
  } else { 
    warning( "length of object != length of op replacement" )
    
  }
  
  x
  
}

#' @rdname formula.parts
#' @aliases op<-,expression-method
#' @name op<-
setReplaceMethod( 'op', 'expression' , .replace.op.plural )


#' @rdname formula.parts
#' @aliases op<-,list-method
#' @name op<-
setReplaceMethod( 'op', 'list' , .replace.op.plural ) 

