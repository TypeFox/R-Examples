simulateFromPrior <-
function (par.prior, times, PRIOR.MODEL=c("none", "OU", "BM", "BMdrift"))
{
	PRIOR.MODEL=match.arg (PRIOR.MODEL)
	
	############################
	# simulate scalar parameters 
	############################
	sdWICluster = runif (1, min=0, max=par.prior$uWICluster)
	sdTSampling = runif (1, min=0, max=par.prior$uTSampling)
	sdResidual = runif (1, min=0, max=par.prior$uResidual)
	
	############################
	# simulate cluster mean 
	############################
	if (PRIOR.MODEL=="none")
	{
		# If not specifying a model, use input mean
		cluster.mean = par.prior$mean
	}
	else 
	{
		# Simulate from a stochastic process
		meanT1 = rnorm (1, mean=par.prior$meanMT1, sd=par.prior$sdMT1)
		sdT1 = runif (1, min=0, max=par.prior$uSDT1)
		meanProc = rnorm (1, mean=par.prior$meanMTProc, sd=par.prior$sdMTProc)
		sdProc = runif (1, min=0, max=par.prior$uSDProc)
		if (PRIOR.MODEL=="OU")
		betaProc = rgamma (1, shape=par.prior$shapeBetaProc, rate=par.prior$rateBetaProc)
		else if (PRIOR.MODEL=="BM")	# Brownian motion without drift
		{
			betaProc = 0
			meanProc = 0
		}
		else	# Brownian motion with drift
		{
			betaProc = 0
		}
		
		cluster.mean = simulateStocProc (diff(times), meanT1=meanT1, sdT1=sdT1, meanProc=meanProc, sdProc=sdProc, betaProc=betaProc, MODEL=PRIOR.MODEL)
	}
	
	return (c (cluster.mean, sdWICluster, sdTSampling, sdResidual))	
}

