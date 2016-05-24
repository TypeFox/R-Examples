hankel.matrix <- function( n, x )
{
###
### This function returns an order n Hankel matrix from the values in
### the order m vector x derived by a circular shift of the values
### in this vector to the left.  Hankel matrices are formed when a given sequence
### of output data and a realization of an underlying state-space or
### hidden Markov model is desired.
###
### arguments
### n = a positive integer value for the order of the Hankel matrix
### x = an order 2 * n + 1 vector of numeric values
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
    if ( m < n )
        stop( "length of argument x is less than n" )
    y <- x
    H <- matrix( 0, nrow=n, ncol=n )
    H[1,] <- y[1:n]
    for ( i in 2:n ) {
        y <- c( y[2:m], y[1] )
        H[i,] <- y[1:n]
    }
    return( H )
}

    