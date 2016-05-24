jacobi.p.polynomials <- function( n, alpha, beta, normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k Jacobi polynomials Pk(a,b,x)
### for orders k=0,1,...n
###
### Parameters
### n = integer highest polynomial order
### alpha = first polynomial parameter
### beta = second polynomial parameter
### normalized = boolean value.  if true, the polynomials are normalized
###
    
    if ( ( abs( alpha ) < 1e-6 ) & ( abs( beta ) < 1e-6 ) )
        return( legendre.polynomials( n, normalized ) )
    if ( abs( alpha - beta ) < 1e-6 ) {
        alpha.prime <- alpha + 0.5
        return( gegenbauer.polynomials( n, alpha.prime, normalized ) )
    }    
    recurrences <- jacobi.p.recurrences( n, alpha, beta, normalized )
    if ( normalized ) {
        ap1 <- alpha + 1
        bp1 <- beta + 1
        abp1 <- alpha + beta + 1
        abp2 <- alpha + beta + 2
        h.0 <- ( 2 ^ abp1 ) * gamma( ap1 ) * gamma( bp1 ) / gamma( abp2 )
        p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
        polynomials <- orthonormal.polynomials( recurrences, p.0 )
    }
    else
        polynomials <- orthogonal.polynomials( recurrences )
    return( polynomials )
}
