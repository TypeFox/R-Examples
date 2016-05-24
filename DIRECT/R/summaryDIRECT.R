summaryDIRECT <-
function (data.name, PERM.ADJUST=FALSE)
{
	file.mcmc.cs = paste (data.name, "_mcmc_cs.out", sep="")
	file.mcmc.pars = paste (data.name, "_mcmc_pars.out", sep="")
	file.mcmc.probs = paste (data.name, "_mcmc_probs.out", sep="")
	file.mcmc.perms = paste (data.name, "_mcmc_perms.out", sep="")

	pars.mcmc = as.matrix (read.delim (file.mcmc.pars, header=FALSE, sep=""))
	cs = as.matrix (read.delim (file.mcmc.cs, header=FALSE, sep=""))
	probs = as.matrix (read.delim (file.mcmc.probs, header=FALSE, sep=""))
	perms = as.matrix (read.delim (file.mcmc.perms, header=FALSE, sep=""))
	if (PERM.ADJUST)
		perms = perms+1
	
	niter = nrow (cs)
	nitem = nrow (probs) / niter
	ncomp = ncol (perms)
	
	alpha.mcmc = cs[,ncol(cs)]
	cs.mcmc = cs[,-ncol(cs)]
	cs.max = apply (cs.mcmc, 1, max)

	probs.perm = array (0, dim=c(niter, nitem, ncomp))
	for (i in 1:niter)
		probs.perm[i,,] = probs[nitem*(i-1)+1:nitem,perms[i,]]
	
	##########################################################
	# Estimate (mean) posterior allocation probability matrix
	##########################################################
	probs.avg = apply (probs.perm, c(2,3), mean)
	colnames (probs.avg) = paste ("C",1:ncol (probs.avg), sep="")
	
	probs.avg.orders = t(apply (probs.avg, 1, order, decreasing=TRUE))
	probs.avg.ordered = t(apply (probs.avg, 1, sort, decreasing=TRUE))
	
	##########################################################
	# Find top two most likely allocations for each item
	##########################################################
	allocprobs = cbind (probs.avg.orders[,1:2], round(probs.avg.ordered[,1:2], digits=3))
	colnames (allocprobs) = c("first", "second", "prob1", "prob2")
	allocprobs = as.data.frame (allocprobs)
	
	##########################################################
	# Assign clusters based on most likely allocations
	##########################################################
	clusters.summary = summary (as.factor (allocprobs$first))
	nclust = nlevels (as.factor (allocprobs$first))
	nclust.all = max (cs.mcmc)
	
	##########################################################
	# Mean posterior estimates of parameters
	##########################################################
	clust.labels = sort(unique (allocprobs$first))
	npars = ncol(pars.mcmc)-2
	pars.perm = array (0, dim=c(niter, nclust.all, npars))
	pars.mcmc.clust = pars.perm
	for (i in 1:niter)
	{
		ind.tmp = which (pars.mcmc[,1]==i)
		pars.tmp = matrix (0, nrow=nclust.all, ncol=npars)
		pars.tmp[1:length(ind.tmp),] = pars.mcmc[ind.tmp, 3:ncol(pars.mcmc)]
		pars.mcmc.clust[i,,] = pars.tmp
		pars.perm[i,,] = pars.tmp[perms[i,],]
	}
	
	pars.avg.mean = apply (pars.perm, c(2,3), mean)
	pars.avg.med = apply (pars.perm, c(2,3), median)
	
	return (list (nitem=nitem, nclust=nclust, top.clust.alloc = allocprobs$first, cluster.sizes=clusters.summary, top.clust.labels = clust.labels, top2allocations=allocprobs, post.alloc.probs = probs.avg, post.clust.pars.mean = pars.avg.mean[clust.labels,], post.clust.pars.median = pars.avg.med[clust.labels,], misc=list (post.pars.mean = pars.avg.mean, post.pars.median = pars.avg.med)))
}

