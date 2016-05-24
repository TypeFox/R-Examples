# -----------------------------------------------------------------------------
# lhs
#   extract and manipulate the left-hand side of R objects.
# -----------------------------------------------------------------------------

#' @include parts.R
NULL


# -------------------------------------
# SINGULAR
# -------------------------------------

#' @rdname formula.parts
#' @aliases .lhs.singular

.lhs.singular <- 
  function(x) 
    if( is.two.sided(x) ) x[[2]] else 
      if( is.one.sided(x) ) NULL   else
        warning( "Could not extract lhs of ", x ) 

 

#' @rdname formula.parts
#' @docType methods
#' @export lhs

setGeneric( 'lhs', function(x, ...) standardGeneric( 'lhs' ) )


# @rdname formula.parts
# @aliases lhs
 
#' @rdname formula.parts
#' @aliases lhs,call-method
setMethod( 'lhs', 'call', .lhs.singular ) 


#' @rdname formula.parts
#' @aliases lhs,formula-method
setMethod( 'lhs', 'formula', .lhs.singular )  

# Note: This is not a replacement method, but rather a method that
# dispatches on the non-standard class '<-'
# @usage \S4method{lhs}{`<-`}(x)

#' @rdname formula.parts
#' @aliases lhs,<--method
setMethod( 'lhs', '<-', function(x) x[[2]] )



# -------------------------------------
# PLURAL
#   Since the 
# -------------------------------------
# setMethod(  'lhs', 'expression', function(x, ... ) lapply( x, lhs, ... ) )

#' @rdname formula.parts
#' @aliases lhs,expression-method
setMethod(  'lhs', 'expression', 
  function(x, ... ) {
    ret <- vector( "expression", length(x) )
    for( i in 1:length(x) ) { 
      lh <- lhs( x[[i]] ) 
      if( ! is.null(lh) )  ret[[i]] <- lh   
    } 
    ret
  }
)

#' @rdname formula.parts
#' @aliases lhs,list-method
setMethod(  'lhs', 'list', function(x, ...) lapply( x, lhs, ... ) )
