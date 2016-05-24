spherical.quadrature.rules <- function( n, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...n.
### An order k quadrature data frame contains the roots and
### abscissa values for the spherical polynomial of degree k
###
### Parameters
### n = integer highest order
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    return( legendre.quadrature.rules( n, normalized ) )
}
