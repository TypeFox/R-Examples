dtrunc <- function( x, spec, a = -Inf, b= Inf, ... )
{
###
### this function computes the density function defined by the spec argument
### for the vector of quantile values in x.  The random variable is truncated
### to be in the interval ( a, b )
###
### Arguments
### x = a numeric vector of quantiles
### spec = a character value for the name of the distribution (e.g., "norm")
### ... = other arguments passed to the corresponding density function
###
    if ( a >= b ) 
        stop( "argument a is greater than or equal to b" )
    tt <- rep( 0, length( x ) )
    g <- get( paste( "d", spec, sep="" ), mode="function" )
    G <- get( paste( "p", spec, sep="" ), mode="function" )
    tt[x >= a & x <= b] <- g( x[x >= a & x <= b], ...) / ( G( b, ... ) - G( a, ... ) )
    return( tt )
}
