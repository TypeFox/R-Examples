ghermite.h.weight <- function( x, mu )
{
###
### This function returns the value of the weight function for
### the generalized Hermite polynomial, Hk(mu,x)
###
### Parameters
### x = a numeric vector function argument
### mu = the polynomial parameter
###
    if ( !is.vector( x ) )
        stop( "function argument is not a vector" )
    if ( mu <= -0.5 )
        stop( "parameter mu is less than or equal to -0.5" )
    if ( mu == 0 )
        return( hermite.h.weight( x ) )
    n <- length( x )
    w <- rep( 0, n )
    mu.times.2 <- 2 * mu
    for ( i in 1:n ) {
        abs.x.i <- abs( x[i] )
        x.i.sq <- x[i] * x[i]
        w[i] <- ( abs.x.i ^ mu.times.2 ) * exp( - x.i.sq )
    }
    return( w )
}
