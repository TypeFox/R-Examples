K.matrix <- function( r, c=r )
{
###
### This function constructs the rc by rc commutation matrix
###
### Arguments
### r = a positive integer value for the number of rows in the H matrices
### c = a positive integer value for the number of columns in the H matrices
###
    if ( missing( r ) )
        stop( "argument r is missing" )
    if ( !is.numeric( r ) )
        stop( "argument r is not numeric" )
    if ( r != trunc( r ) )
        stop( "argument r is not an integer" )
    if ( r < 2 )
        stop( "argument r is less than 2" )
    if ( !is.numeric( c ) )
        stop( "argument c is not numeric" )
    if ( c != trunc( c ) )
        stop( "argument c is not an integer" )
    if ( c < 2 )
        stop( "argument c is less than 2" )
    H <- H.matrices( r, c )
    p <- r * c
    K <- matrix( 0, nrow=p, ncol=p )
    for ( i in 1:r ) {
        for ( j in 1:c ) {
            Hij <- H[[i]][[j]]
            K <- K + ( Hij %x% t( Hij ) )
        }
    }
    return( K )
}
