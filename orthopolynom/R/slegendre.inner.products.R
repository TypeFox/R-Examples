slegendre.inner.products <- function( n )
{
###
###	This function returns a vector with n+1 elements
###	containing the inner product of an order k shifted Legendre polynomial,
###	Pstar-k(x), with itself( i.e. the norm squared) for orders k=0,1,...,n
###
###	Parameter
###	n = integer highest polynomial order
###
	if ( n < 0 )
		stop( "negative highest polynomial order" )
	if ( n != round( n ) )
		stop( "highest polynomial order is not integer" )
	inner.products <- rep( 0, n + 1 )
	j <- 1
	for ( k in 0:n ) {
		inner.products[j] <- 1 / ( 2 * k + 1 )
		j <- j + 1
	}
	return ( inner.products )
}	
