ptrunc <- function( q, spec, a = -Inf, b = Inf, ... )
{
###
### this function computes the distribution function defined by the spec argument
### for the vector of quantile values in x.  The random variable is truncated
### to be in the interval ( a, b )
###
### Arguments
### q = a numeric vector of quantiles
### spec = a character value for the name of the distribution (e.g., "norm")
### ... = other arguments passed to the corresponding density function
###
    if ( a >= b )
        stop( "argument a is greater than or equal to b" )
    tt <- q
    aa <- rep( a, length( q ) )
    bb <- rep( b, length( q ) )
    G <- get( paste( "p", spec, sep="" ), mode="function" )
    tt <- G( apply( cbind( apply( cbind( q, bb ), 1, min ), aa ), 1, max ), ... )
    tt <- tt - G ( aa, ... )
    result <- tt / ( G( bb, ... ) - G ( aa, ... ) )
    return( result )
}

