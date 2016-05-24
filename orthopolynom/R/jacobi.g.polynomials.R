jacobi.g.polynomials <- function( n, p, q, normalized=FALSE )
{
###
###   This function returns a list with n+1 elements
###   containing the order k Jacobi polynomials Gk(p,q,x)
###   for orders k=0,1,...,n
###
###   Parameters
###   n = integer highest polynomial order
###   p = first polynomial parameter
###   q = second polynomial parameter
###   normalized = boolean value.  if true, the polynomials are normalized
###
    recurrences <- jacobi.g.recurrences( n, p, q, normalized )
    if ( normalized ) {
        h.0 <- gamma( q ) * gamma( p - q + 1 ) / gamma( p + 1 )
        p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
        polynomials <- orthonormal.polynomials( recurrences, p.0 )
    }
    else
        polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
