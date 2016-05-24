matrix.inverse <- function( x )
{
###
### this function returns the inverse of a square matrix
###
### Parameters
### x = a square numeric matrix
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    return( solve( x ) )
}
