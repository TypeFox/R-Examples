computeLogLikInd <-
function (cluster.par, gene.data, n.time, n.rep)
{
	cluster.mean = cluster.par[1:n.time]
	sdWICluster = cluster.par[n.time+1]
	sdTSampling = cluster.par[n.time+2]
	sdResidual = cluster.par[n.time+3]
	
# compute mean vector
	mean.vec = rep (cluster.mean, times=n.rep)
	
# compute covariance matrix
	Smatrix = matrix (1, nrow=n.time, ncol=n.time)
	Imatrix = matrix (0, nrow=n.time, ncol=n.time)
	diag (Imatrix) = 1
	diagonal.block = sdWICluster^2 * Smatrix + (sdTSampling^2 + sdResidual^2) * Imatrix
	offdiagonal.block = sdWICluster^2 * Smatrix + sdTSampling^2 * Imatrix
	
	cov.matrix = matrix (0, nrow=n.time*n.rep, ncol=n.time*n.rep)
	for (r in 1:n.rep)
	{
		cov.matrix[(n.time * (r-1) + 1:n.time), (n.time * (r-1) + 1:n.time)] = diagonal.block 
	}
	for (r1 in 1:(n.rep-1))
	{
		for (r2 in (r1+1):n.rep)
		{
			cov.matrix[(n.time * (r1-1) + 1:n.time), (n.time * (r2-1) + 1:n.time)] = offdiagonal.block 
			cov.matrix[(n.time * (r2-1) + 1:n.time), (n.time * (r1-1) + 1:n.time)] = offdiagonal.block 
		}
	}
	
# compute likelihoods
	loglik = dMVNorm (gene.data, mean=mean.vec, sigma=cov.matrix, log=TRUE)
	
	return (loglik)
}

