chebyshev.t.weight <- function( x )
{
###
### This function returns the value of the weight function
### for the Chebyshev polynomial of the first kind, Tk( x )
###
### Parameter
### x = the function argument
###
    n <- length( x )
    y <- rep( 0, n )
    for ( i in 1:n ) {
        if ( ( x[i] > -1 ) && ( x[i] < +1 ) )
            y[i] <- 1 / ( sqrt( 1 - x[i] * x[i] ) )
    }
    return( y )
}
