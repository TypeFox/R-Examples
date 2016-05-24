chebyshev.s.weight <- function( x )
{
###
### This function returns the value of the weight function
### for the Chebyshev polynomial of the second kind, Sk( x )
###
### Parameter
### x = the function argument
###
    n <- length( x )
    y <- rep( 0, n )
    for ( i in 1:n ) {
        if ( ( x[i] > -2 ) && ( x[i] < 2 ) )
            y[i] <- ( sqrt( 1 - 0.25 * x[i] * x[i] ) )
    }
    return( y )
}
