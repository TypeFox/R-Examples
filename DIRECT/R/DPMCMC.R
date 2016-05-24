DPMCMC <-
function (file.mcmc.cs, file.mcmc.pars, file.mcmc.probs, file.size, data, SKIP, nTime, times, c.curr, par.prior, PAR.INIT=FALSE, sdWICluster.curr=0.5, sdTSampling.curr=0.5, sdResidual.curr=0.5, alpha.curr, alpha.prior.shape, alpha.prior.rate, sd.prop, nIter, burn.in, step.size, nRepeat=1, nResample, seed.value, RNORM.METHOD=c("chol", "eigen", "svd"), SAMPLE.C=c("FRBT", "Neal"), PRIOR.MODEL=c("none", "OU", "BM", "BMdrift"), ALPHA.METHOD=c("Gibbs", "MH"), OUTPUT.CLUST.SIZE=FALSE, PRINT=FALSE)
{
	set.seed (seed.value)
	
	RNORM.METHOD = match.arg (RNORM.METHOD)
	SAMPLE.C = match.arg (SAMPLE.C)
	PRIOR.MODEL = match.arg (PRIOR.MODEL)
#	PROPOSAL = match.arg (PROPOSAL)
	ALPHA.METHOD = match.arg (ALPHA.METHOD)    # Gibbs or MH update for concentration parameter alpha
	
	# re-organize data into array of nGene by nTime by nRep
	nGene = nrow (data)
	nRep = (ncol (data) - SKIP) / nTime
	nPar = nTime + 3
	
	ts = array (0, dim = c(nGene, nTime, nRep))
	for (r in 1:nRep)
	{
		ts[,,r] = as.matrix (data[,SKIP + (0:(nTime-1))*nRep + r])
	}
	
	ts.mtx = as.matrix (data[,(SKIP+1): ncol (data)])

	####################################
	# initialization
	####################################
	need.append = 0
		
	c.curr = reassignLabelsVec (c.curr)
	nClust = nlevels (as.factor (c.curr))
	# generate initial values of cluster-specific parameter values
	cluster.mean = computeClusterMean (ts, c.curr, nClust)	# result is a matrix of C by nTime
	if (PAR.INIT)
	{
		sdWICluster.curr = runif (nClust, min=0, max=par.prior$uWICluster)
		sdTSampling.curr = runif (nClust, min=0, max=par.prior$uTSampling)
		sdResidual.curr = runif (nClust, min=0, max=par.prior$uResidual)
	}
	
	par.unique.curr = matrix (0, nrow=nClust, ncol=nPar)
	for (l in 1:nClust)
		par.unique.curr[l,] = c (as.vector (cluster.mean[l,]), sdWICluster.curr[l], sdTSampling.curr[l], sdResidual.curr[l])
	
	need.update = 1
	####################################
	# MH updates
	####################################
	for (k in 1:nIter)
	{
#		cat ("iteration:", k, "\n")
 		cat (k, "")
		if (SAMPLE.C=="FRBT")
			result = updateClusterMemberships_MH_FRBT (c.curr=c.curr, par.unique.curr=par.unique.curr, alpha.curr=alpha.curr, data=ts.mtx, nGene=nGene, nTime=nTime, nRep=nRep, par.prior=par.prior, times=times, PRIOR.MODEL=PRIOR.MODEL, PRINT=PRINT)
		else
			result = updateClusterMemberships_MH_Neal (c.curr=c.curr, par.unique.curr=par.unique.curr, alpha.curr=alpha.curr, data=ts.mtx, nGene=nGene, nTime=nTime, nRep=nRep, nRepeat=nRepeat, par.prior=par.prior, times=times, PRIOR.MODEL=PRIOR.MODEL, PRINT=PRINT)
		
		result = reassignLabels (result$c.curr, result$par.unique.curr)
		c.curr = result$c
		par.unique.curr = result$pars
		c.counts = table (as.factor (c.curr))
		nClust = length (c.counts)

#		cat ("update alpha\n")
		########################################
		# update Dirichlet parameter alpha
		# now Gibbs updating
		########################################
		alpha.new = updateDAlpha (alpha.curr, alpha.prior.shape, alpha.prior.rate, sd.prop$alpha, nClust, nGene, ALPHA.METHOD)
#		cat ("alpha:", alpha.new, "\n")

#		cat ("done updating alpha\n")
		########################################
		# update parameter value of each cluster
		########################################
#		cat ("updating cluster-specific parameter values\n")
		par.unique.new = matrix (0, nrow=nClust, ncol=nPar)

		# update mean profile of each gene
#		cat ("updating mean profile of each gene\n")
		mean.ind = matrix (0, nrow=nGene, ncol=nTime)
		for (i in 1:nGene)
		{
			par.curr = par.unique.curr[c.curr[i],]
			mean.ind[i,] = updateIndMean (ts[i,,], par.curr[1:nTime], par.curr[nTime+1], par.curr[nTime+2], par.curr[nTime+3], RNORM.METHOD=RNORM.METHOD)
		}
			
		
#		cat ("done updating mean profile of each gene\n")
		cluster.size.vec = rep ("NA", nGene)
		# for each unique c
#		cat ("updating parameters\n")
		for (c in 1:nClust)
		{
#			cat (c)
#			cat ("cluster:", c, "\n")
			# update cluster mean profile
			gene.id = which (c.curr==c)
			cluster.size = length (gene.id)
#			cat (paste ("(",cluster.size,")", sep=""))
			cluster.size.vec[c] = cluster.size

			cluster.mean = par.unique.curr[c, 1:nTime]
			sdWICluster = par.unique.curr[c, nTime+1]
			sdTSampling = par.unique.curr[c, nTime+2]
			sdResidual = par.unique.curr[c, nTime+3]

			if (cluster.size>1)
				cluster.mean.update = as.vector (apply (mean.ind[gene.id,], 2, mean))
			else
				cluster.mean.update = as.vector(mean.ind[gene.id,])
#			cat ("done cluster means\n")

			# update standard deviation terms
			sdRET.curr = par.unique.curr[c,nTime+1:3]
			for (j in 1:3)
			{
				if (cluster.size==1 & j==1)
					sdRET.curr[1] = 0
				else
				{
					sdRET.update = updateRandomEffectsInd (i=j, sdRET.curr=sdRET.curr, ts.cluster=ts.mtx[gene.id,], cluster.mean=cluster.mean.update, n.time=nTime, n.rep=nRep, upper=par.prior[[j]], sdProp=sd.prop[[j]])
					sdRET.curr = sdRET.update
				}
			}
#			cat ("done sd residual\n")

			par.unique.new[c, c(1:(nTime+3))] = c (cluster.mean.update, sdRET.update)
		}
#		cat ("\n")
		
		par.unique.curr = par.unique.new
		alpha.curr = alpha.new
		
		##########################################################
		# output labels and cluster-specific parameter estimates
		##########################################################
		if (OUTPUT.CLUST.SIZE)
		{
			if (k==1)
				write (cluster.size.vec, file.size, ncolumns=nGene, sep="\t")
			else
				write (cluster.size.vec, file.size, ncolumns=nGene, sep="\t", append=TRUE)
		}
		
		result.cs = as.vector (c(c.curr,alpha.curr))
		if (PRINT)
			print (result.cs)
		ncol.output = length (result.cs)
		if ((k > (nIter * burn.in)) && ((k %% step.size) == 0))
		{
			if (need.append == 0)
			{
				count = 1
				write (result.cs, file.mcmc.cs, ncolumns = ncol.output, sep = "\t")
				result.pars = cbind (rep(count, nrow (par.unique.curr)), 1:nrow(par.unique.curr), par.unique.curr)
				write.table (result.pars, file.mcmc.pars, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t")
				
				count = count + 1
				need.append <- 1
			}
			else
			{
				write (result.cs, file.mcmc.cs, ncolumns = ncol.output, append = TRUE, sep = "\t")
				result.pars = cbind (rep(count, nrow (par.unique.curr)), 1:nrow(par.unique.curr), par.unique.curr)
				write.table (result.pars, file.mcmc.pars, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, sep = "\t")
				
				count = count + 1
			}
		}
	}

	cat ("\nMCMC step finished\n")
}

