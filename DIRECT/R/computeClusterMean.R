computeClusterMean <-
function (ts, c.curr, nclust)
{
	c.counts = summary (as.factor (c.curr), maxsum=length (c.curr))
	c.unique = unique (c.curr)
#	c.unique.sorted = c.unique[order(c.unique, decreasing=FALSE)]
	result = matrix (0, nrow=nclust, ncol=dim(ts)[2])
	for (l in 1:nclust)
	{
		if (c.counts[l]>1)
			result[l,] = as.vector (apply (ts[which (c.curr==c.unique[l]),,], 2, mean))
		else 
			result[l,] = as.vector (apply (ts[which (c.curr==c.unique[l]),,], 1, mean))
	}

	return (result)
}

