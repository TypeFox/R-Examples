updateIndMean <-
function (data.i, cluster.mean, sdWICluster, sdTSampling, sdResidual, RNORM.METHOD)
{
	d = nrow (data.i)
	n = ncol (data.i)
	Smatrix = matrix (1, nrow=d, ncol=d)
	Imatrix = Smatrix - 1
	diag (Imatrix) = 1

	prior.precision = solve (sdWICluster^2 * Smatrix + sdTSampling^2 * Imatrix) # are there analytic results?
	prior.mean = cluster.mean

	sample.precision = sdResidual^(-2) * Imatrix
	sample.mean = apply (data.i, 1, mean)

	results = rmvnormpost (n, prior.mean, prior.precision, sample.mean, sample.precision, RNORM.METHOD)
	return (results)
}

