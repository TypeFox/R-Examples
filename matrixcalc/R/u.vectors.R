u.vectors <- function( n )
{
###
### This function constructs an identity matrix I of order
### p = n * ( n +1 ) / 2.  It also builds a lower triangular square matrix of
### order n.  The value of element [i,j] is the column number
### in the identify matrix.  It is the mapping from the coordinates
### to the column vector in the identity matrix.
###
### argument
### n = a positive integer value
###
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    p <- n * ( n + 1 ) / 2
    I <- diag( rep( 1, p ) )
    k <- matrix( 0, nrow=n, ncol=n )
    for ( j in 1:n ) {
        for ( i in j:n ) {
            k[i,j] <- ( j - 1 ) * n + i - 0.5 * j * ( j - 1 )
        }
    }
    return( list( k=k, I=I ) )
}
