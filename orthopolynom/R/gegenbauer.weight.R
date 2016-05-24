gegenbauer.weight <- function( x, alpha )
{
###
### This function returns the value of the weight function
### for the Gegenbauer polynomial, Ck( alpha, x )
###
### Parameters
### x = the function argument
### alpha = the polynomial parameter
###
    n <- length( x )
    y <- rep( 0, n )
    pow <- alpha - 0.5
    for ( i in 1:n ) {
        if ( ( x[i] > -1 ) && ( x[i] < 1 ) ) {
            y[i] <- ( 1 - x[i] * x[i] ) ^ pow
        }
    }
    return( y )
}
