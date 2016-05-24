schebyshev.t.polynomials <- function( n, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k shifted Chebyshev polynomials of the first kind, Tk(x),
### for orders k=0,1,...n
###
### Parameters
### n = integer highest polynomial order
### normalized = boolean value.  if true, the polynomials are normalized
###
    require( polynom )
    recurrences <- schebyshev.t.recurrences( n, normalized )
    if ( normalized ) {
        h.0 <- pi
        p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
        polynomials <- orthonormal.polynomials( recurrences, p.0 )
    }
    else
        polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
