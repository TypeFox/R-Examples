gegenbauer.inner.products <- function( n, alpha )
{
###
### This function returns a vector with n+1 elements
### containing the inner product of an order k Gegenbauer polynomial, Ck(alpha,x),
### with itself (i.e. norm squared) for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### alpha = polynomial parameter
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
        return( legendre.inner.products( n ) )
###
### alpha = 1.0
### special case is the Chebyshev polynomial of the second kind U
###
    if ( abs( alpha - 1.0 ) < 1e-6 )
        return( chebyshev.u.inner.products( n ) )
###
### initialize vector for the inner products
###
    inner.products <- rep( 1, n + 1 )
###
### alpha = 0.0
### special case is related to the Chebyshev polynomial of the first kind T
###
    if ( abs( alpha ) < 1e-6 ) {
        inner.products[1] <- pi
        j <- 2
        for ( k in 1:n ) {
            inner.products[j] <- ( 2 * pi ) / ( k ^ 2 )
            j <- j + 1
        }
        return( inner.products )
    }
###
### general case
###
    coef <- pi * ( 2 ^ ( 1 - 2 * alpha ) )
    j <- 1
    for ( k in 0:n ) {
        num <- coef * gamma( k + 2 * alpha )
        den <- factorial( k ) * ( k + alpha ) * ( gamma( alpha ) ) ^ 2
        inner.products[j] <- num / den
        j <- j + 1
    }
    return( inner.products )
}
