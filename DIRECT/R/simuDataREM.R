simuDataREM <-
function (pars.mtx, dt, T, ntime, nrep, nsize, times, method=c("eigen", "svd", "chol"), model=c("OU", "BM"))
{
	method = match.arg (method)
	model = match.arg (model)
	
	data.tmp = simuMultiClustREM (pars.mtx=pars.mtx, dt=dt, T=T, ntime=ntime, nrep=nrep, ngene=nsize, times=times, method=method, model=model)
	
	labels = sample (1:nrow(data.tmp$data), nrow (data.tmp$data), replace=FALSE)
	data.scrambled = data.tmp$data[labels,]
	
	return (list (means=data.tmp$means, data=data.scrambled))
}

