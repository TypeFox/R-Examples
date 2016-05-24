.loadSgnU <- function (x)
{														##	change the signs of the columns of a loadings matrix such, that each column's absolute maximum is positive
	idx.max <- apply (abs (x), 2, which.max)
	sgn <- sign (x[cbind (idx.max, 1:ncol (x))])
	if (length (sgn) == 1)
		return (x * sgn)
	return (x %*% diag (sgn))
}
