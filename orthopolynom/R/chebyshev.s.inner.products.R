chebyshev.s.inner.products <- function( n )
{
###
###	This function returns a vector with n+1 elements
###	containing the inner product of an order k Chebyshev polynomial
###	of the second kind, Sk(x) with itself (i.e. the norm squared)
###	for orders k=0,1,...n
###
###	Parameter
###	n = integer highest polynomial order
###
	if ( n < 0 )
		stop( "negative highest polynomial order" )
	if ( n != round( n ) )
		stop( "highest polynomial order is not integer" )
	inner.products <- rep( pi, n + 1 )
	return ( inner.products )
}
