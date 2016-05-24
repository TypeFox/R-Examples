chebyshev.c.quadrature.rules <- function( n, normalized = FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...,n.
### An order k quadrature data frame contains the roots and
### abscissa values for the order k Chebyshev polynomial
### of the first kind, Ck(x), and of degree k
###
### Parameter
### n = integer highest order
### normalized = boolean value.  If TRUE, recurrences are for normalized polynomials
###
    if ( n <= 0 )
        stop( "highest order is not positive" )
    if ( n != round( n ) )
        stop( "highest order is not an integer" )
    recurrences <- chebyshev.c.recurrences( n, normalized )
    inner.products <- chebyshev.c.inner.products( n )
    return( quadrature.rules( recurrences, inner.products ) )
}   
