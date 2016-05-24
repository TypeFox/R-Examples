legendre.quadrature.rules <- function( n, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...n.
### An order k quadrature data frame contains the roots
### abscissas values for the Legendre polynomial of degree k
###
### Parameters
### n = integer highest order
### normalized = a boolean value.  if true, the recurrences are for normalized polynomials
###
    recurrences <- legendre.recurrences( n, normalized )
    inner.products <- legendre.inner.products( n )
    return( quadrature.rules( recurrences, inner.products ) )
}
