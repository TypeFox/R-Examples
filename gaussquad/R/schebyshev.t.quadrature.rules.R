schebyshev.t.quadrature.rules <- function( n, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...,n.
### An order k quadrature data frame contains the roots and
### abscissa values for the order k shifted Chebyshev polynomial
### of the first kind, Tk(x), and of degree k
###
### Parameter
### n = integer highest order
### normalized = boolean value.  If TRUE, recurrences are for normalized polynomials
###
    recurrences <- schebyshev.t.recurrences( n, normalized )
    inner.products <- schebyshev.t.inner.products( n )
    return( quadrature.rules( recurrences, inner.products ) )
}   
