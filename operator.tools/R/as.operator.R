# ---------------------------------------------------------------------
# as.operator:
#   coerce to an operator
#   similar to 'as.function', but will only work for operators.
# ---------------------------------------------------------------------

as.operator <- function(x, ...) {
  UseMethod( 'as.operator', x ) 
}

as.operator.function <- function(x, ...)
  if( is.operator(x) ) return(x) else
    stop( x, " cannot be coerced to an operator." )


as.operator.character <- function(x, ...)
  if( x %in% operators(...) ) eval( as.name(x) ) else
    stop( x, " cannot be coerced to an operator." )


as.operator.name <- function(x, ... ) 
  if( deparse(x) %in% operators(...) ) eval( x ) else
    stop( x, " cannot be coerced to an operator." )  


