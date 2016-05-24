creation.matrix <- function( n )
{
###
### this function returns an order n creation matrix that is useful
### in numerical mathematics such as the solution of differential equations.
###
### argument
### n = a positive integer greater than 1
###
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( !is.vector( n ) )
        stop( "argument n is not the proper data type" )
    if ( length( n ) > 1 )
        stop( "argument n is not a scalr number" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    H <- shift.down( diag( seq( 1, n ) ) )
    return( H )
}
