upper.triangle <- function( x )
{
###
### this function returns the lower triangular matrix portion of matrix x
###
### Parameters
### x = a numeric matrix
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square numeric matrix" )
    y <- x
    y[row(x) > col(y)] <- 0
    return( y )
}
