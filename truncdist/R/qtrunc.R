qtrunc <- function( p, spec, a = -Inf, b = Inf, ... )
{
###
### This function evaluates the inverse of the truncated distribution function
### or quantile function for the given vector of probabilities and
### the specified distribution
###
### Arguments
### p = a numeric vector of probabilities
### spec = a character value for the name of the distribution (e.g., "norm")
### ... = other arguments passed to the corresponding density function
###
    if ( a >= b )
        stop( "argument a is greater than or equal to b" )
    tt <- p
    G   <- get( paste( "p", spec, sep="" ), mode = "function" )
    Gin <- get( paste( "q", spec, sep="" ), mode = "function" )
    result <- pmin( pmax( a, Gin( G( a, ... ) + 
              p * ( G( b, ... ) - G( a, ... )), ... ) ), b )
    return( result )
}    
    
