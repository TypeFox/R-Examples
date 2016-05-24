ddirichlet <- function (x, alpha) 
{
    if (length(x) != length(alpha))
	stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
        logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
        s <- sum((alpha - 1) * log(x))
        pd <- exp(sum(s) - logD)
	  pd[any(x < 0 | x > 1)] <- 0
	  if(sum(x) != 1) pd <- 0
	  return(pd)
}