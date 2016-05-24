monic.polynomials <- function( monic.recurrences )
{
###
###	This function returns a list with n+1 elements
###	containing the order k monic polynomials for orders k=0,1,...,n
###
###	Parameter
###	monic.recurrences = a data frame containing the parameters a and b
###
	require( polynom )
	np1 <- nrow( monic.recurrences )
	n <- np1 - 1
	polynomials <- as.list( rep( NULL, np1 ) )
	p.0 <- polynomial( c(1) )
	polynomials[[1]] <- p.0
	a <- monic.recurrences$a
	b <- monic.recurrences$b
	j <- 0
	while( j < n ) {
		aj <- a[j+1]
		bj <- b[j+1]
		monomial <- polynomial( c( -aj, 1 ) )
		if ( j == 0 ) {
			p.jp1 <- monomial
		}
		else {
			p.jm1 <- polynomials[[j]]
			p.j   <- polynomials[[j+1]]
			p.jp1 <- monomial * p.j - bj * p.jm1
		}
		polynomials[[j+2]] <- p.jp1
		j <- j + 1
	}
	return( polynomials )
}
