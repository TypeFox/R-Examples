spherical.recurrences <- function( n, normalized=FALSE )
{
###
###	This function returns a data frame with n+1 rows and four columns
###	containing the coefficients c, d, e and f of the recurrence relations
###	for the order k spherical polynomial, Pk(x), and orders k=0,1,...,n
###
###	Parameter
###	n = integer highest polynomial order
###	normalized = a boolean value.  if true, the recurrences are for normalized polynomials
###
	return( legendre.recurrences( n, normalized ) )
}	
