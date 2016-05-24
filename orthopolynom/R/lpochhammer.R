lpochhammer <- function( z, n )
{
###
###	This function returns the logarithm of Pochhammer's symbol
###
###	Parameters
###	z = the argument of the symbol
###	n = integer number of terms in the symbol
###
	if ( n < 0 )
		stop( "n is negative" )
	else if ( n == 0 )
		return ( 1 )
	else
		return ( lgamma( z + n ) - lgamma( z ) )
	return ( NULL )
}
