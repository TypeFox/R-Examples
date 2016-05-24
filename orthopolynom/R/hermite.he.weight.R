hermite.he.weight <- function( x )
{
###
### This function returns the value of the weight function
### for the scaled Hermite polynomial, He-k( x )
###
### Parameter
### x = function argument
###
    return( exp( -0.5 * ( x * x ) ) )
}
