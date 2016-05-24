jacobi.p.weight <- function( x, alpha, beta )
{
###
### This function returns the value of the weight function
### for the Jacobi polynomial, Pk( alpha, beta, x )
###
### Parameters
### x = the function argument
### alpha = the first polynomial parameter
### beta = the second polynomial parameter
###
    if ( alpha < -1 )
        stop( "argument alpha is less than -1" )
    if ( beta < -1 )
        stop( "argument beta is less than =1" )
    if ( ( abs( alpha ) < 1e-6 ) & ( abs( beta ) < 1e-6 ) )
        return( legendre.weight( x ) )
    n <- length( x )
    y <- rep( 0, n )
    for ( i in 1:n ) {
        if ( ( x[i] > -1 ) && ( x[i] < +1 ) ) {
            y[i] <- ( ( 1 - x[i] ) ^ alpha ) * ( ( 1 + x[i] ) ^ beta )
        }    
    }
    return ( y )
}
