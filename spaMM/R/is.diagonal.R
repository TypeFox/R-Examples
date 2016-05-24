# is.square.diagonal <- function( x, tol=1e-8 ) {
#   if (inherits(x,"ddiMatrix")) return(isDiagonal( x )) ## implies square in Matrix package
#   if ( ! is.matrix( x ) ) stop( "argument x is not a numeric matrix" )
#   if ( ncol(x) != nrow(x) ) return( FALSE )
#   if ( ! is.numeric( x ) ) stop( "argument x is not a numeric matrix" )
#   y <- x
#   diag( y ) <- rep( 0, nrow( y ) )
#   return( max(abs( y )) < tol )
# }
## isDiagonal must work on matrix as well as Matrix


is.identity <- function( x, matrixcheck=FALSE, tol=1e-8 ) {
  if (inherits(x,"Matrix")) {
    if (inherits(x,"ddiMatrix") ) return(x@diag=="U")
    if (! isDiagonal( x ) ) return( FALSE )
    ## hence now diagonal:
    return(max(abs(range(diag(x))-1))<tol)
  } else {
    if (matrixcheck) {
      return(ncol(x)==nrow(x) && max(abs(x- diag(ncol(x))))<tol)
    } else return(FALSE)  
  }
}

