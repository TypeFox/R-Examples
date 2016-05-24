jacobi.p.quadrature.rules <- function( n, alpha, beta, normalized = FALSE )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...,n
### An order k quadrature data frame contains the roots and
### abscissa values for the Jacobi polynomial, Pk(alpha,beta,x), and
### of degree k
###
### Parameters
### n = integer highest order
### alpha = first polynomial parameter
### beta = second polynomial parameter
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    recurrences <- jacobi.p.recurrences( n, alpha, beta, normalized )
    inner.products <- jacobi.p.inner.products( n, alpha, beta )
    return( quadrature.rules( recurrences, inner.products ) )
}   
