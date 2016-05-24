jacobi.g.quadrature.rules <- function( n, p, q, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...,n
### An order k quadrature data frame contains the roots and
### abscissa values for the Jacobi polynomial, Gk(a,b,x), and
### of degree k
###
### Parameters
### n = integer highest order
### p = first polynomial parameter
### q = second polynomial parameter
###
    recurrences <- jacobi.g.recurrences( n, p, q, normalized )
    inner.products <- jacobi.g.inner.products( n, p, q )
    return( quadrature.rules( recurrences, inner.products ) )
}
