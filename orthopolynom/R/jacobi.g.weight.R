jacobi.g.weight <- function( x, p, q )
{
###
###   This function returns the value of the weight function
###   for the Jacobi polynomial, Gk( p, q, x )
###
###   Parameters
###   x = the function argument
###   p = the first polynomial parameter
###   q = the second polynomial parameter
###
    n <- length( x )
    y <- rep( 0, n )
    for ( i in 1:n ) {
       if ( ( x[i] > 0 ) && ( x[i] < 1 ) ) {
          t1 <- ( 1 - x[i] ) ^ ( p - q )
          t2 <- ( x[i] ) ^ ( q - 1 )
          y[i] <- t1 * t2
       }
    }
    return( y )
}
