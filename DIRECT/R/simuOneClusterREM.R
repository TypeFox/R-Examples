simuOneClusterREM <-
function (pars, dt, T, ntime, nrep, ngene, times, method=c("eigen", "svd", "chol"), model=c("OU", "BM"))
{
	method = match.arg (method)
	model = match.arg (model)
	
	start = rnorm (1, mean=pars[4], sd=pars[5])
	ts.tmp = simuProcess (start=start, beta=pars[8], mu=pars[6], sigma=pars[7], dt=dt, T=T, model=model)
	
	sdWIcluster = pars[1]
	sdTSampling = pars[2]
	sdResidual = pars[3]
	
	covmtx = computeCovMtx (sdWIcluster, sdTSampling, sdResidual, ntime, nrep)
	data.tmp = rMVNorm (ngene, mean=rep (ts.tmp[times+1, 2], each=nrep), sigma=covmtx, method=method)
	
	return (list (mean=ts.tmp[times+1,2], data=data.tmp))
	
}

