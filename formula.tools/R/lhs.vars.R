# -----------------------------------------------------------------------------
# lhs.vars
# -----------------------------------------------------------------------------

#' @aliases lhs.vars
#' @rdname get.vars
#' @export

setGeneric( 'lhs.vars', function(x, ... ) standardGeneric( 'lhs.vars' ) )

#' @rdname get.vars
#' @aliases .lhs.vars 
.lhs.vars <- function(x, ..., data=NULL) 
{
  if( 
      class( x[[1]] )   == 'name' &&
      deparse( x[[1]] ) %in% operators() 
  ) {
    get.vars( lhs(x), ..., data=data ) 
  } else {
    warning( "There is no relational operator defined for ", deparse(x)  )
  }

}


#' @rdname get.vars
#' @aliases lhs.vars,formula-method
setMethod( 'lhs.vars' , 'formula', .lhs.vars )

#' @rdname get.vars
#' @aliases lhs.vars,call-method
setMethod( 'lhs.vars' , 'call'   , .lhs.vars )

#' @rdname get.vars
#' @aliases lhs.vars,expression-method
setMethod( 'lhs.vars' , 'expression', function(x,...) lapply(x, .lhs.vars, ...))



