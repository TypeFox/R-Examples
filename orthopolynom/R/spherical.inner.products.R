spherical.inner.products <- function( n )
{
###
###	This function returns a vector with n+1 elements
###	containing the inner product of an order k Spherical polynomial, Pk(x),
###	with itself (i.e. the norm squared) for orders k=0,1,...n
###
###	Parameter
###	n = integer highest polynomial order
###
	return ( legendre.inner.products( n ) )
}	
