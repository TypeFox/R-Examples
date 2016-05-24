ultraspherical.quadrature.rules <- function( n, alpha, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...,n.
### An order k quadrature data frame contains the roots and
### abscissa values for the ultraspherical polynomial of degree k
###
### Parameters
### n = integer highest order
### alpha = polynomial parameter
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    return( gegenbauer.quadrature.rules( n, alpha, normalized ) )
}	
