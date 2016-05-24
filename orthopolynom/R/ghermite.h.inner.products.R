ghermite.h.inner.products <- function( n, mu )
{
###
### This function returns a vector with n+1 elements
### containing the inner product of an order k generalized Hermite polynomial, Hk(mu,x),
### with itself (i.e. the norm squared) for orders k=0,1,...,n
###
### Parameter
### n = integer highest polynomial order
### mu = first parameter
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( mu <= -0.5 )
        stop( "mu less than or equal to -0.5" )
    if ( abs( mu ) < 1e-6 )
        return( hermite.h.inner.products( n ) )
    inner.products <- rep( 1, n + 1 )
    sqrt.pi <- sqrt( pi )
    log.2 <- log( 2 )
    mu.plus <- mu + 0.5
    j <- 1
    for ( k in 0:n ) {
        floor.k <- floor( k / 2 )
        floor.kp1 <- floor( ( k + 1 ) / 2 )
        log.inner.product <- (2 * k ) * log.2 + lgamma( floor.k + 1 ) + lgamma( floor.kp1 + mu.plus )
        inner.products[j] <- exp( log.inner.product )
        j <- j + 1
    }   
    return ( inner.products )
}

