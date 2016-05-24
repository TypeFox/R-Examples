###################################################
#		All Permutation indexes
###################################################
Permu.Ind <- function(y, npermu, pair=F)
{
	res <- NULL
	if (!pair)	## unpaired data
	{
		len <- length(y)
		if (factorial(len) > npermu)	## random samples
		{
			res <- matrix(0, npermu, len)
			for (i in 1 : npermu)
			{
				res[i, ] <- sample(c(1:len), len, replace=F)
			}
		} else	## all samples
		{
			res <- matrix(0, factorial(len), len)
			res1 <- permn(1 : len)
			for (i in 1 : factorial(len))
			{
				res[i, ] <- res1[[i]]
			}
		}
	} else	## paired data
	{
		len <- length(y) / 2
		if (2^len > npermu)	## random samples
		{
			res <- matrix(0, npermu, len*2)
			for (i in 1 : npermu)
			{
				res1 <- sample(c(0, 1), len, replace=T)
				res[i, (1:len)*2-1] <- (1:len)*2 - res1
				res[i, (1:len)*2] <- (1:len)*2 - (1-res1)
			}
		} else	## all samples
		{
			res <- matrix(0, 2^len, len*2)
			res1 <- hcube(rep(2, len)) - 1
			for (i in 1 : (2^len))
			{
				res[i, (1:len)*2-1] <- (1:len)*2 - res1[i,]
				res[i, (1:len)*2] <- (1:len)*2 - (1-res1[i,])
			}
		}
	}
	
	return(res)
}
