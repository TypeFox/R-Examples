vech <- function( x )
{
###
### returns a stack of the lower triangular matrix as a matrix with 1 column
###
### Parameters
### x = a numeric matrix square matrix
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square numeric matrix" )
    return( t( t( x[!upper.tri(x)] ) ) )
}
