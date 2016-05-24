pochhammer <- function( z, n )
{
###
###	This function returns the value of Pochhammer's symbol
###	calculated as ( z ) * ( z + 1 ) * ... * ( z + n - 1 )
###
###	Parameters
###	z = the argument of the symbol
###	n = integer number of terms in the product
###
	if ( n < 0 )
		stop( "n is negative" )
	else if ( n == 0 )
		return ( 1 )
	else {
		result <- 1
		for ( i in 1:n ) {
			result <- result * ( z + i - 1 )
		}
		return ( result )
	}
	return ( NULL )
}
