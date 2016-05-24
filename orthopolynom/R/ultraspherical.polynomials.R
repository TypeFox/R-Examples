ultraspherical.polynomials <- function( n , alpha, normalized=FALSE )
{
###
###	This function returns a list with n+1 elements
###	containing the order k ultraspherical polynomials, Ck(alpha,x),
###	for orders k=0,1,...,n
###
###	Parameters
###	n = integer highest polynomial order
###
	require( polynom )
	polynomials <- gegenbauer.polynomials( n, alpha, normalized )
	return( polynomials )
}
