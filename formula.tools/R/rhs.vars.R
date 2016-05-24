# -----------------------------------------------------------------------------
# rhs.vars
# -----------------------------------------------------------------------------

#' @rdname get.vars
#' @aliases rhs.vars
#' @export
setGeneric( 'rhs.vars', function(x, ... ) standardGeneric( 'rhs.vars' ) )

#' @rdname get.vars
#' @aliases .rhs.vars
.rhs.vars <- 
  function(x, ..., data=NULL) 
  {
    if( 
        class( x[[1]] )   == 'name' &&
        deparse( x[[1]] ) %in% operators()  
    ) {
  
      term.rhs <- terms( x, data=data, ... ) 
      labels   <- attr( term.rhs, 'term.labels' )
      order    <- attr( term.rhs, 'order' )
      vars.rhs <- labels[ order == 1 ]
  
      vars.rhs
  
    } else {
      warning( "There is no relational operator defined for ", deparse(x)  )
    }
  
  }





#' @rdname get.vars
#' @aliases rhs.vars,formula-method
setMethod( 'rhs.vars' , 'formula', .rhs.vars )

#' @rdname get.vars
#' @aliases rhs.vars,call-method
setMethod( 'rhs.vars' , 'call'   , .rhs.vars )

#' @rdname get.vars
#' @aliases rhs.vars,expression-method
setMethod( 'rhs.vars' , 'expression', function(x,...) lapply(x, rhs.vars, ...))



