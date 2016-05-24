laguerre.inner.products <- function( n )
{
###
###	This function returns a vector with n+1 elements
###	containing the inner product of order k Laguerre polynomials, Lk(x),
###	with itself (i.e. the norm squared) for orders k=0,1,...,n
###
###	Parameter
###	n = integer highest polynomial order
###
	return( glaguerre.inner.products( n, 0 ) )
}
