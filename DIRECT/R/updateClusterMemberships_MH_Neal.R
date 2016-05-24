updateClusterMemberships_MH_Neal <-
function (c.curr, par.unique.curr, alpha.curr, data, nGene, nTime, nRep, nRepeat, par.prior, times, PRIOR.MODEL=c("none", "OU", "BM", "BMdrift"), PRINT=FALSE)
{
	PRIOR.MODEL = match.arg (PRIOR.MODEL)
	nPar = ncol (par.unique.curr)
	need.update = 1
	for (i in 1:nGene)
	{
#		cat (i)
		
		#######################################
		# update c_i R times
		#######################################
		for (r in 1:nRepeat)
		{
			if (need.update)
			{
				# relabel the cluster memberships;
				# and reorder the cluster-specific parameter estimates
				# so rows of the matrix are in the same order as after reassignment
				result = reassignLabels (c.curr, par.unique.curr)
				c.curr = result$c
				par.unique.curr = result$pars
				c.counts = table (as.factor (c.curr))
				nClust = length (c.counts)
			}
			
			# check 0 in par.unique.curr
			par.row.sum = apply (par.unique.curr, 1, sum)
			
			if (length (which (par.row.sum==0)))
			{
				print ("after relabeling")
				#cat ("k=", k, "\n")
				cat ("i=", i, "\n")
				print (dim (par.unique.curr))
				cat ("rows", which (par.row.sum==0), "\n")
				stop ("zeros!")
			}
			
			###########################################
			# generate a proposal c*_i from Eq. (5.4)
			###########################################
			c.counts.noi = c.counts
			c.counts.noi[c.curr[i]] = c.counts.noi[c.curr[i]] - 1
			prob.i = c (c.counts.noi / (nGene-1+alpha.curr), alpha.curr / (nGene-1+alpha.curr))
			cdf.i = diffinv (prob.i)[-1]
			ci.new = simulateObserved (cdf.i)
			
			c.new = c.curr
			
			#############################################
			# if c*_i is new,
			# sample new value from prior distribution,
			# calculate acceptance probability
			# and decided whether to accept c*_i
			#############################################
			# if c*_i is new
			if (ci.new > nClust)
				par.i.new = simulateFromPrior (par.prior, times, PRIOR.MODEL)
			else
				par.i.new = par.unique.curr[ci.new, ]
			
			par.row.sum = apply (par.unique.curr, 1, sum)
			
			if (length (which (par.row.sum==0)))
			{
				print ("before relabeling")
				#cat ("k=", k, "\n")
				cat ("i=", i, "\n")
				print (dim (par.unique.curr))
				cat ("rows", which (par.row.sum==0), "\n")
				stop ("zeros!")
			}
			
			##################################################
			# compute acceptance probability as in Eq. (5.3)
			##################################################
			clusterParHastingsRatio = computeClusterParHastingsRatio (par.i.new, par.unique.curr[c.curr[i],], as.vector (data[i,]), nTime, nRep)
			accept.prob.log = min (0, clusterParHastingsRatio)
			
			########################################
			# decide whether to accept c*_i or not
			########################################
			if (log (runif (1)) < accept.prob.log)
			{
				c.new[i] = ci.new
				if (ci.new > nClust)
					par.unique.curr = rbind (par.unique.curr, par.i.new)
				need.update = 1
			}
			else
				need.update = 0
			
			if (PRINT) cat ("r=", r, ":", c(c.new, par.i.new), "\n")
			
			c.curr = c.new
		}
		
	}
	if (PRINT) cat ("\n")
	
	return (list (c.curr=c.curr, par.unique.curr=par.unique.curr))
}

