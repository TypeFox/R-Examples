normalizeProbMtx <-
function (mtx, direc)
{
	margin.sum = as.vector (apply (mtx, direc, sum))
	if (direc==1)
	{
		result = mtx / margin.sum
	}
	if (direc==2)
	{
		result = t(t(mtx)/margin.sum)
	}
	return (result)
}

