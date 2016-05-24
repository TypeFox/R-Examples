laguerre.polynomials <- function( n , normalized=FALSE )
{
###
### This function returns a list with n+1 elements
### containing the order k Laguerre polynomials, Lk(x),
### for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
###
    polynomials <- glaguerre.polynomials( n, 0, normalized )
    return( polynomials )
}
