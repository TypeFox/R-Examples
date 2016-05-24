computeLogLik <-
function (cluster.par, ts.cluster, n.time, n.rep, DATA.FORMAT=0)
{
    if (is.vector(ts.cluster)) {
        ts.cluster <- matrix(ts.cluster, ncol = length(ts.cluster))
    }
	
	cluster.mean = cluster.par[1:n.time]
	sdWICluster = cluster.par[n.time+1]
	sdTSampling = cluster.par[n.time+2]
	sdResidual = cluster.par[n.time+3]
	
	# compute mean vector
	if (DATA.FORMAT==0)
		mean.vec = rep (cluster.mean, each=n.rep)
	else
		mean.vec = rep (cluster.mean, n.rep)
	
	# compute covariance matrix
	cov.matrix = computeCovMtx (sdWICluster, sdTSampling, sdResidual, n.time, n.rep, DATA.FORMAT)
	
	# compute likelihoods
	loglik = sum(dMVNorm (ts.cluster, mean=mean.vec, sigma=cov.matrix, log=TRUE))
	
	return (loglik)
}

