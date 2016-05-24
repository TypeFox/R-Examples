rtrunc <- function( n, spec, a = -Inf, b = Inf, ... )
{
###
### This function generates n random numbers that are drawn from
### the specified truncated distribution.
###
### Arguments
### n = a positive integer for the number of random numbers generated
### spec = a character value for the name of the distribution (e.g., "norm")
### ... = other arguments passed to the corresponding density function
###
    if ( a >= b )
        stop( "argument a is greater than or equal to b" )
    u <- runif( n, min = 0, max = 1 )
    x <- qtrunc( u, spec, a = a, b = b, ... )
    return( x )
}
