extrunc <- function( spec, a = -Inf, b = Inf, ... )
{
###
### This function computes the expected value of a truncated random variable
### with the given specification over the truncated domain using
### numerical integration
###
### Arguments
### spec = a character value to specify the distribution
### a = a numeric value for the lower bound of the truncation interval
### b = a numeric value for the upper bound of the truncation interval
### ... = one or more arguments used by the density function
###
    if ( a >= b )
        stop( "argument a is greater than or equal to b" )
    f <- function( x ) { x * dtrunc( x, spec, a = a, b = b, ... ) }
    return( integrate( f, lower = a, upper = b )$value )
}
