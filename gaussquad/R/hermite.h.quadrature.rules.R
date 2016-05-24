hermite.h.quadrature.rules <- function( n, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...,n
### An order k quadrature data frame contains the roots and
### abscissa values for the degree k normalized Hermite polynomial
###
### Parameter
### n = integer highest order
### normalized = a boolean value.  If true, recurrences are for normalized polynomials
###
    recurrences <- hermite.h.recurrences( n, normalized )
    inner.products <- hermite.h.inner.products( n )
    return( quadrature.rules( recurrences, inner.products ) )
}
