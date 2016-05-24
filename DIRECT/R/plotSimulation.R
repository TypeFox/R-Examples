plotSimulation <-
function (simudata, times=1:ntime, nsize, ntime=length(times), nrep, skip=0, ...)
{
	nclust = length (nsize)
	for (k in 1:nclust)
	{
		simudata.ts = arrayTSData (simudata$data[which (simudata$data[,1]==k),-1], N=nsize[k], T=ntime, R=nrep, skip=skip)
		
		simudata.med = apply (simudata.ts, c(1,2), mean)
		matplot (times, t(simudata.med), ...)
		lines (times, simudata$means[k,-1], col="red")
	}
}

