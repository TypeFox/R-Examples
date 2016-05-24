rDirichlet <- function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

dDirichlet <- function (x, alpha, log=FALSE)
{
	dlog = lgamma (sum (alpha)) + sum ((alpha-1)*log(x)) - sum (lgamma (alpha))
	result = ifelse (!log, exp (dlog), dlog)
	return (result)
}

