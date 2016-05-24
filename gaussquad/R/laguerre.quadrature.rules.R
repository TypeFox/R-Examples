laguerre.quadrature.rules <- function ( n, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...n.
### An order k quadrature data frame contains the roots and
### abscissa values for the Laguerre polynomial of degree k.
###
### Parameter
### n = integer highest order
### normalized = a boolean value.  If true, recurrences are for normalized polynomials
###
    return( glaguerre.quadrature.rules( n, 0, normalized ) )
}
