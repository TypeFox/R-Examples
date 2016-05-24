glaguerre.inner.products <- function( n, alpha )
{
###
###	This function returns a vector with n+1 elements
###	containing the inner product of an order k generalized Laguerre polynomial,
###	Lk(alpha,x ) with itself (i.e. the norm squared) for orders k=0,1,...n
###
###	Parameter
###	n = integer highest polynomial order
###	alpha = polynomial parameter
###
	if ( n < 0 )
		stop( "negative highest polynomial order" )
	if ( n != round( n ) )
		stop( "highest polynomial order is not integer" )
	if ( alpha <= -1 )
		stop( "alpha less than or equal to -1" )
	inner.products <- rep( 0, n + 1 )
	j = 1;
	for ( k in 0:n ) {
		inner.products[j] <- exp( lgamma( alpha + k + 1 ) - lfactorial( k ) )
		j <- j + 1
	}	
	return ( inner.products )
}
