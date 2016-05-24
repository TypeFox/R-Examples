ultraspherical.recurrences <- function( n, alpha, normalized=FALSE )
{
###
###	This function returns a data frame with n+1 and four columns
###	containing the coefficients c, d, e and f of the recurrence relations
###	for the order k ultraspherical polynomials, Ck(alpha,x),
###	and for orders k=0,1,...,n
###
###	Parameters
###	n = integer highest polynomial order
###	normalized = a boolean value.  If true, recurrences are for normalized polynomials
###
	return( gegenbauer.recurrences( n, alpha, normalized ) )
}
