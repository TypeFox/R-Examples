glaguerre.weight <- function( x, alpha )
{
###
###	This function returns the value of the weight function
###	for the generalized Laguerre polynomial, Lk( x, alpha )
###
###	Parameters
###	x = function argument
###	alpha = polynomial argument
###
	n <- length( x )
	y <- rep( 0, n )
	for ( i in 1:n ) {
		if ( x[i] > 0 )
			y[i] <- exp(-x[i]) * (x[i] ^ alpha )
	}
	return( y )
}
