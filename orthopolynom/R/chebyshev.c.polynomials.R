chebyshev.c.polynomials <- function( n, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k Chebyshev polynomials of the first kind, Ck(x),
### for orders k=0,1,...n
###
### Parameters
### n = integer highest polynomial order
### normalized = boolean value.  if true, the polynomials are normalized
###
    recurrences <- chebyshev.c.recurrences( n, normalized )
    if ( normalized ) {
        h.0 <- 8 * pi
        p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
        polynomials <- orthonormal.polynomials( recurrences, p.0 )
    }
    else
        polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
