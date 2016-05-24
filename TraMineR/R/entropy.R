## Compute the entropy of a distribution

entropy <- function(distrib, base=exp(1))
	{
		distrib <- distrib[distrib!=0]
		p <- distrib/sum(distrib)
		e <- -sum(p*log(p, base=base))
		return(e)
	}
