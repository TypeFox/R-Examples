jacobi.p.inner.products <- function( n, alpha, beta )
{
###
### This function returns a vector with n+1 elements
### containing the inner product of an order k Jacobi polynomial
### Pk(alpha,beta,x) for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### a = first parameter
### b = second parameter
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( alpha <= -1 )
        stop( "alpha less than or equal to -1" )
    if ( beta <= -1 )
        stop( "beta less than or equal to -1" )
    if ( ( abs( alpha ) < 1e-6 ) & ( abs( beta ) < 1e-6 ) )
        return( legendre.inner.products( n ) )
    if ( abs( alpha - beta ) < 1e-6 ) {
        alpha.prime <- alpha + 0.5
        return( gegenbauer.inner.products( n, alpha.prime ) )
    }    
    ab <- alpha + beta
    abp1 <- alpha + beta + 1
    ap1 <- alpha + 1
    bp1 <- beta + 1
    coef <- 2 ^ abp1
    inner.products <- rep( 1, n + 1 )
    j <- 1
    for ( k in 0:n ) {
        num <- coef * gamma( k + ap1 ) * gamma( k + bp1 )
        den <- ( 2 * k + abp1 ) * factorial( k ) * gamma( k + abp1 )
        inner.products[j] <- num / den
        j <- j + 1
    }
    return( inner.products )
}
