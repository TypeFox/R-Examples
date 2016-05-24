glaguerre.polynomials <- function( n, alpha, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k generalized Laguerre polynomial Lk(alpha,x),
### for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### alpha = polynomial parameter
### normalized = boolean value.  if true, the polynomials are normalized
###
    recurrences <- glaguerre.recurrences( n, alpha, normalized )
    if ( normalized ) {
        h.0 <- gamma( alpha + 1 )
        p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
        polynomials <- orthonormal.polynomials( recurrences, p.0 )
    }
    else
        polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
