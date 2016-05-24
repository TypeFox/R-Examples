hermite.h.weight <- function( x )
{
###
### This function returns the value of the weight function
### for the Hermite polynomial, Hk( x )
###
### Parameter
### x = function argument
###
    return( exp( - ( x * x ) ) )
}
