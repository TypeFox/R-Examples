gegenbauer.quadrature.rules <- function( n, alpha, normalized=FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...,n
### An order k quadrature data frame contains the roots and
### abscissa values for the Gegenbauer polynomial of degree k
###
### Parameters
### n = integer highest order
### alpha = polynomial parameter
### normalized = boolean value.  If TRUE, recurrences are for normalized polynomials
###
### alpha = 0.5
### special case is the Legendre polynomial
###
    if ( abs( alpha - 0.5 ) < 1e-6 )
        return( legendre.quadrature.rules( n, normalized ) )
###
### alpha = 1.0
### special case is the Chebyshev polynomial of the second kind U
###
    if ( abs( alpha - 1.0 ) < 1e-6 )
        return( chebyshev.u.quadrature.rules( n, normalized ) )
###
### all other cases including alpha = 0
###
    recurrences <- gegenbauer.recurrences( n, alpha, normalized )
    inner.products <- gegenbauer.inner.products( n, alpha )
    return( quadrature.rules( recurrences, inner.products ) )
}
