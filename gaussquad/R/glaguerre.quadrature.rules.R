glaguerre.quadrature.rules <- function( n, alpha, normalized=FALSE )
{
###
###	This function returns a list with n elements
###	containing the order k quadrature rule data frames
###	for orders k=1,2,...,n.
###	An order k quadrature data frame contains the roots and
###	abscissa values for the generalized Laguerre polynomial of degree k.
###
###	Parameters
###	n = integer highest order
###	alpha = polynomial parameter
###     normalized = a boolean value.  if true, the recurrences are for normalized polynomials
###
	recurrences <- glaguerre.recurrences( n, alpha, normalized )
	inner.products <- glaguerre.inner.products( n, alpha )
	return( quadrature.rules( recurrences, inner.products ) )
}
