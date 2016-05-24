toeplitz.matrix <- function( n, x )
{
###
### this function returns the order n Toeplitz matrix derived from
### the order 2 * n - 1 vector
###
### arguments
### n = a positive integer value for the order of the Hankel matrix
### x = an order 2 * n - 1 vector of numeric values
###
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( n != trunc( n ) )
        stop( "argument n ix not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    if ( !is.vector( x ) )
        stop( "argument x is not a vector" )
    m <- length( x )
    if ( m != ( 2 * n - 1 ))
        stop( "length of argument x is not comparable with n n" )
    T <- matrix( 0, nrow=n, ncol=n )
    for ( i in 1:n ) {
        for ( j in 1:n ) {
            T[i,j] <- x[(i - j) + n]
        }
    }
    return( T )
}
