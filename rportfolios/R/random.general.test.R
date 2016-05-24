random.general.test <- function( n=2, k=n, p=0.5, x.u=1/k )
{
###
### This function generates a general random portfolio of n investments.
### There are k non-zero long and short positions.
### The probability that a position is a long position is p.  The maximum absolute
### exposure for any investment is x.u.  The default value is 1 / k.
###
### Arguments
### n = a positive integer for the number of investments in the portfolio
### p = a positive numeric value for the probability that an investment has positive weight
### x.u = a positive numeric value for the maximum absolute exposure to an investment
###
    if ( n <= 0 ) {
        stop( "Argument n is not positive" )
    }    
    if ( ( p < 0 ) || ( p > 1 ) ) {
        stop( "Argument p is not a valid probability" )
    }
    if ( x.u <= 0 ) {
        stop( "Argument x.u is not positive" )
    }
    if ( k < n ) {
        answer <- random.general.test( n=k, k, p, x.u )
        weight <- answer$x
        investments <- sample( 1:n, k, replace=FALSE )
        x <- rep( 0, n )
        x[investments] <- weight
        result <- list( x=x, iter=1 )
        return( result )
    }
    signs <- 2 * rbinom(n, 1, p ) - 1
    exposures <- x.u * runif( n )
    x <- signs * exposures
    result <- list( x=x, iter=1 )
    return( result )
}
