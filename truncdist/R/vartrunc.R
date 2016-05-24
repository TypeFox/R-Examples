vartrunc <- function( spec, a = -Inf, b = Inf, ... )
{
###
### This function computes the variance of the truncated random variable
###
### Arguments
### spec = a character value to specify the distribution
### a = a numeric value for the lower bound of the truncation interval
### b = a numeric value for the upper bound of the truncation interval
### ... = one or more arguments used by the density function
###
    if ( a >= b )
        stop( "argument a is greater than or equal to b" )
    ex <- extrunc( spec, a = a, b = b, ... )
    f <- function( x ) { ( x - ex )^2 * dtrunc( x, spec, a = a, b = b, ... ) }
    tt <- integrate( f, lower = a, upper = b )$value 
    return( tt )
}
