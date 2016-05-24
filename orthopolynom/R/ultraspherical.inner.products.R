ultraspherical.inner.products <- function( n, alpha )
{
###
###	This function returns a vector with n+1 elements
###	containing the inner product of an order k ultraspherical polynomial
###	with itself (i.e. norm squared) for orders k=0,1,...,n
###
###	Parameters
###	n = integer highest polynomial order
###	alpha = polynomial parameter
###
	return( gegenbauer.inner.products( n, alpha ) )
}	
