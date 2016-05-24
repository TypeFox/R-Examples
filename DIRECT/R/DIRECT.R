DIRECT <-
function (data, data.name="Output", SKIP=0, nTime, times=1:nTime, c.curr, uWICluster=1, uTSampling=1, uResidual=1, meanVec=rep(0, nTime), meanMT1=0, sdMT1=0.2, meanMTProc=0, sdMTProc=0.5, uSDT1=0.2, uSDProc=1, shapeBetaProc=0.5, rateBetaProc=0.5, PAR.INIT=TRUE, sdWICluster.curr=0.5, sdTSampling.curr=0.5, sdResidual.curr=0.5, alpha.curr=0.01, alpha.prior.shape=0.01, alpha.prior.rate=1, WICluster.prop.sd = 0.2, TSampling.prop.sd = 0.2, Residual.prop.sd = 0.2, alpha.prop.sd=0.2, nIter, burn.in, step.size, nRepeat=1, nResample, seed.value, RNORM.METHOD=c("chol", "eigen", "svd"), SAMPLE.C=c("FRBT", "Neal"), PRIOR.MODEL=c("none", "OU", "BM", "BMdrift"), ALPHA.METHOD=c("Gibbs", "MH"), RELABEL.THRESHOLD=0.01, OUTPUT.CLUST.SIZE=FALSE, PRINT=FALSE)
{
	##############################################
	# Initialization
	##############################################
	RNORM.METHOD = match.arg (RNORM.METHOD)
	SAMPLE.C = match.arg (SAMPLE.C)
	PRIOR.MODEL = match.arg (PRIOR.MODEL)
#	PROPOSAL = match.arg (PROPOSAL)
	ALPHA.METHOD = match.arg (ALPHA.METHOD)    # Gibbs or MH update for concentration parameter alpha

	par.prior = list (uWICluster=uWICluster, uTSampling=uTSampling, uResidual=1, mean=meanVec, meanMT1=meanMT1, sdMT1=sdMT1, meanMTProc=meanMTProc, sdMTProc=sdMTProc, uSDT1=uSDT1, uSDProc=uSDProc, shapeBetaProc=shapeBetaProc, rateBetaProc=rateBetaProc)

	sd.prop = list (WICluster = WICluster.prop.sd, TSampling = TSampling.prop.sd, Residual = Residual.prop.sd, alpha=alpha.prop.sd)

	#######################################################
	# re-organize data into array of nItem by nTime by nRep
	#######################################################
	nItem = nrow (data)
	nRep = (ncol (data) - SKIP) / nTime
	
	ts = array (0, dim = c(nItem, nTime, nRep))
	for (r in 1:nRep)
	{
		ts[,,r] = as.matrix (data[,SKIP + (0:(nTime-1))*nRep + r])
	}
	
	#####################################################################
	# Markov Chain Monte Carlo to sample partitions and parameters 
	# under the Dirichlet-process prior
	#####################################################################
	file.mcmc.cs = paste (data.name, "_mcmc_cs.out", sep="")
	file.mcmc.pars = paste (data.name, "_mcmc_pars.out", sep="")
	file.mcmc.probs = paste (data.name, "_mcmc_probs.out", sep="")
	file.mcmc.perms = paste (data.name, "_mcmc_perms.out", sep="")
	file.size = paste (data.name, "_mcmc_size.out", sep="")
	cat ("start MCMC...\n")
	DPMCMC (file.mcmc.cs=file.mcmc.cs, file.mcmc.pars=file.mcmc.pars, file.mcmc.probs=file.mcmc.probs, file.size=file.size, data=data, SKIP=SKIP, nTime=nTime, times=times, c.curr=c.curr, par.prior=par.prior, PAR.INIT=TRUE, alpha.curr=alpha.curr, alpha.prior.shape=alpha.prior.shape, alpha.prior.rate=alpha.prior.rate, sd.prop=sd.prop, nIter=nIter, burn.in=burn.in, step.size=step.size, nRepeat=nRepeat, nResample=nResample, seed.value=seed.value, RNORM.METHOD=RNORM.METHOD, SAMPLE.C=SAMPLE.C, PRIOR.MODEL=PRIOR.MODEL, ALPHA.METHOD=ALPHA.METHOD, OUTPUT.CLUST.SIZE=OUTPUT.CLUST.SIZE, PRINT=PRINT)

	#####################################################################
	# Resample MCMC realizations to estimate 
	# posterior allocation probability matrix
	#####################################################################
	# read in mcmc samples
	pars.mcmc = as.matrix (read.delim (file.mcmc.pars, header=FALSE))
	cs = as.matrix (read.delim (file.mcmc.cs, header=FALSE))
	cs.mcmc = cs[, -ncol(cs)]
	alpha.mcmc = cs[,ncol(cs)]
	max.c = max (cs.mcmc)
	
	cat ("maximum number clusters across stored MCMC samples:", max.c, "\n")
	cat ("start resampling of mixing proportions...\n")
	resampleClusterProb (file.mcmc.probs, ts, nItem, nTime, nRep, pars.mcmc, cs.mcmc, alpha.mcmc, nstart=1, nres=nResample)
	cat ("\nresampling step finished\n")
	
	##########################################################################
	# Relabel (Algorithm 2 in Stephens 2000 JRSSB) to handle label-switching
	##########################################################################
	cat ("start relabeling step...\n")
	probs.mcmc = as.matrix (read.delim (file.mcmc.probs, header=FALSE))
#	probs = as.double (as.vector (t(probs.mcmc)))
#	nIter = as.integer (nrow (cs.mcmc))
#	result = .C ("relabel_R", perm_mtx = as.integer (rep (0, length=max.c*nIter)), probs, nIter, as.integer(nItem), as.integer(max.c), RELABEL.THRESHOLD, 0)
#	perms.mtx = matrix (result$perm_mtx, byrow=TRUE, ncol=max.c)
	perms.mtx = relabel (probs.mcmc, nIter=nrow(cs.mcmc), nItem=nItem, nClust=max.c, RELABEL.THRESHOLD)
	write.table (perms.mtx, file.mcmc.perms, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	cat ("relabeling step finished\n")
}

