matrix.power <- function( x, k )
{
###
### this function computes the k-th power of order n square matrix x
### if k is zero, the order n identity matrix is returned.  argument k
### must be an integer
###
### arguments
### x = a numeric square matrix
### k = an integer value exponent
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not numeric" )
    if ( !is.numeric( k ) )
        stop( "argument k is not numeric" )
    if ( !is.vector( k ) )
        stop( "argument k is not the proper data type" )
    if ( length( k ) > 1 )
        stop( "argument k is not a scalr number" )
    if ( k != trunc( k ) )
        stop( "argument k is not an integer" )
    n <- nrow( x )
    if ( k == 0 ) {
        return( diag( 1, n ) )
    }    
    if ( k < 0 ) {
        return( matrix.power( matrix.inverse(x), -k ) )
    }
    x.k <- x
    i <- 2
    while ( i <= k ) {
        x.k <- x.k %*% x
        i <- i + 1
    }
    return( x.k )
}
