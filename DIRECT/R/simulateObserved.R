simulateObserved <-
function (cdf)
{
	n = length (cdf)
	tmp = runif (1, min=0, max=1)
	obs = 1
	for (k in 2:n)
	{
		if (tmp>cdf[k-1] & tmp<=cdf[k])
			obs = k
	}
	return (obs)	
}

