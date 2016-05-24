arrayTSData <-
function (data, N, T, R, skip)
{
	ts = array (0, dim = c(N, T, R))
	for (r in 1:R)
		ts[,,r] = as.matrix (data[,skip + (0:(T-1))*R + r])

	return (ts)
}

