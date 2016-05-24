gegenbauer.polynomials <- function( n, alpha, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k Gegenbauer polynomial Ck(alpha,x),
### for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### alpha = polynomial parameter
### normalized = boolean value.  if true, the polynomials are normalized
###
    if ( n < 0 )
        stop( "highest polynomial order is less than zero" )
    if ( n != round(n) )
        stop( "highest polynomial order is not an integer" )
    if ( alpha <= -0.5 )
        stop( "alpha is less than or equal to -0.5" )
###
### alpha = 0.5
### special case is the Legendre polynomial
###
    if ( abs( alpha - 0.5 ) < 1e-6 )
        return( legendre.polynomials( n, normalized ) )
###
### alpha = 1.0
### special case is the Chebyshev polynomial of the second kind U
###
    if ( abs( alpha - 1.0 ) < 1e-6 )
        return( chebyshev.u.polynomials( n, normalized ) )
###
### all other cases including alpha = 0
###
    recurrences <- gegenbauer.recurrences( n, alpha, normalized )
    if ( normalized ) {
        if ( abs( alpha ) < 1e-6 ) {
            h.0 <- pi
        }
        else {
            h.0 <- sqrt( pi ) * gamma( alpha + 0.5 ) / gamma( alpha + 1 )
        }    
        p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
        polynomials <- orthonormal.polynomials( recurrences, p.0 )
    }
    else
        polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
