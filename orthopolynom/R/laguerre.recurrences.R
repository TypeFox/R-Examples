laguerre.recurrences <- function( n, normalized=FALSE )
{
###
###	This function returns a data frame with n+1 and four columns
###	containing the coefficients c, d, e and f of the recurrence relations
###	for the order k Laguerre polynomials, Lk(x),
###	and for orders k=0,1,...,n
###
###	Parameters
###	n = integer highest polynomial order
###	normalized = a boolean value.  If true, recurrences are for normalized polynomials
###
	return( glaguerre.recurrences( n, 0, normalized ) )
}
