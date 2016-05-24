simuMultiClustREM <-
function (pars.mtx, dt, T, ntime, nrep, ngene, times, method=c("eigen", "svd", "chol"), model=c("OU", "BM"))
{
	nclust = nrow (pars.mtx)									# number of clusters
	result = matrix (0, nrow=sum(ngene), ncol=ntime*nrep+1)		# first column cluster membership
	means = matrix (0, nrow=nclust, ncol=ntime+1)				# first column cluster membership
	
	method = match.arg (method)
	model = match.arg (model)

	for (k in 1:nclust)
	{
		print (k)
		data.tmp = simuOneClusterREM (pars=pars.mtx[k,], dt=dt, T=T, ntime=ntime, nrep=nrep, ngene=ngene[k], times=times, method=method, model=model)
		means[k,] = c(k, data.tmp$mean)
		if (k==1)
			result[1:ngene[k],] = cbind (rep (k, ngene[k]), data.tmp$data)
		else
			result[sum(ngene[1:(k-1)])+1:ngene[k],] = cbind (rep (k, ngene[k]), data.tmp$data)
	}

	return (list (means=means, data=result))
}

