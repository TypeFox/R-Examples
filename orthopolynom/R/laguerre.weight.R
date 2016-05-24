laguerre.weight <- function( x )
{
###
###	This function returns the value of the weight function
###	for the Laguerre polynomial, Lk( x )
###
###	Parameter
###	x = function argument
###
	n <- length( x )
	y <- rep( 0, n )
	for ( i in 1:n ) {
		if ( x[i] > 0 )
			y[i] <- exp(-x[i])
	}
	return( y )
}
