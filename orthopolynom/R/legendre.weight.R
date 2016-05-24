legendre.weight <- function( x )
{
###
###	This function returns the value of the weight function
###	for the Legendre polynomial, Pk( x )
###
###	Parameter
###	x = the function argument
###
	n <- length( x )
	y <- rep( 0, n )
	for ( i in 1:n ) {
		if ( ( x[i] > -1 ) && ( x[i] < +1 ) )
			y[i] <- 1
	}
	return( y )
}
