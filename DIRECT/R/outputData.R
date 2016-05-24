outputData <-
function (datafilename, parfilename, meanfilename, simudata, pars, nitem, ntime, nrep)
{
	if (is.vector (pars))
	{
		nclust=1
		pars = matrix (pars, ncol=length(pars))
		truemean = matrix (simudata$mean, ncol=length (simudata$mean))
	}
	else
	{
		nclust=nrow(pars)
		truemean = simudata$mean
	}
	
	# write simulated data to output
	write (t(simudata$data), datafilename, sep="\t", ncolumns=ncol(simudata$data))
	
	# write simulation parameters to output
	write (nitem, parfilename)
	write (ntime, parfilename, append=TRUE)
	write (nrep, parfilename, append=TRUE)
	if (ncol(truemean)>ntime)
		truepars = cbind (truemean, matrix(pars[,1:3],ncol=3))
	else
		truepars = cbind (rep(1:nclust), truemean, matrix(pars[,1:3],ncol=3))
	write.table (truepars, parfilename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	
	# write mean profiles to output
	simudata.ts = arrayTSData (simudata$data, N=nitem, T=ntime, R=nrep, skip=ifelse(nclust==1, 0, 1))
	simudata.mean = apply (simudata.ts, c(1,2), mean)
	write.table (simudata.mean, meanfilename, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

