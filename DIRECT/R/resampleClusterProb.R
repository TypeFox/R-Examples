resampleClusterProb <-
function (file.out, ts, nitem, ntime, nrep, pars.mcmc, cs.mcmc, alpha.mcmc, nstart=1, nres=1)
{
#	set.seed (seed.value)
	
	nclust = max (cs.mcmc)
#	alpha.ncol = ncol (cs.mcmc)
	niter = nrow (cs.mcmc)
	parsncol = ncol(pars.mcmc)
	
	cs.max = apply (cs.mcmc, 1, max)
	
	for (r in nstart:niter)
	{
#		cat ("r=",r,"\n")
		cat (r,"")
		
		result = matrix (0, nrow=nitem, ncol=nclust)
		# No need for resampling if there is only one cluster in the current MCMC sample
		if (cs.max[r] == 1)
		{
			result[,1] = 1
		}
		else
		{
			normden = matrix (-Inf, nrow=nitem, ncol=nclust)
			for (c in 1:nclust)
			{
#				cat ("c=",c,"\n")
				# Find index corresponding to cluster c in sample r
				cind = which (pars.mcmc[,1]==r & pars.mcmc[,2]==c)
				if (length(cind)>0)
				{
					# Extract cluster-specific parameter estimates
					pars = as.vector (pars.mcmc[cind, 3:parsncol])
					# For each gene, compute log likelihood 
					# given current cluster assignment and parameter estimates
					for (i in 1:nitem)
						normden[i,c] = computeLogLikInd (pars, as.vector (ts[i,,]), ntime, nrep)
				}
			}
			
			# Resample nres times for each MCMC sample
			# Only mixing proportions change each time
			for (q in 1:nres)
			{
#				cat ("q=",q,"\n")
				# Generate mixing proportions from Dirichlet prior
				# (prior allocation probabilities)
				# given current number of clusters
#				w = as.vector (rDirichlet (1, alpha=rep(alpha.mcmc[r] / cs.max[r], cs.max[r])))
				w = as.vector (rDirichlet (1, alpha=rep(alpha.mcmc[r], cs.max[r])))
				if (cs.max[r] < nclust)
					w = c (w, rep(0, nclust-cs.max[r]))
				
				# Posterior probability of allocation is proportional to 
				# the product of likelihood and prior probability of allocation		
#				cat ("dim of t(normdne): ", dim (t(normden)), "\n")
#	 			cat ("length of log (w): ", length (log(w)), "\n")
#				print (log(w))
				joint = t(t(normden) + log (w))
				result = result + normalizeProbMtx (exp(joint), 1)
			}
			result = result / nres
		}
		
		# Write posterior allocation matrix to output file
		if (r==nstart)
			write (t(result), file.out, ncolumns=ncol(result), sep="\t")
		else
			write (t(result), file.out, ncolumns=ncol(result), append=TRUE, sep="\t")
	}
}

