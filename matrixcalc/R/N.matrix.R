N.matrix <- function( n )
{
###
### This function returns the n*n by n*n matrix that is the sum of the
### implicit commutation matrix Kn and an identity matrix
###
### Arguments
### n = a positive integer value for the underlying matrix
###
    if ( missing( n ) )
        stop( "argument n is missing" )
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    return( ( K.matrix( n ) + diag( rep( 1, n * n ) ) ) / 2 )
}
