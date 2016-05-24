updateClusterMemberships_MH_FRBT <-
function (c.curr, par.unique.curr, alpha.curr, data, nGene, nTime, nRep, par.prior, times, PRIOR.MODEL=c("none", "OU", "BM", "BMdrift"), PRINT=FALSE)
{
#	PROPOSAL = match.arg (PROPOSAL)
	PRIOR.MODEL = match.arg (PRIOR.MODEL)
	
	nPar = ncol (par.unique.curr)
	need.update = 1
	for (i in 1:nGene)
	{
		if (PRINT) cat (i, ",")
		
		#############################################
		# find size of parameter of current cluster
		#############################################
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
			if (PRINT) print (c.counts)
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
			
		
		#############################################################
		# generate a proposal c*_i uniformly 
		# from 1 through nGene+1, excluding gene i 
		# or from 1 through nClust+1, excluding c_i
		#############################################################
#		if (PROPOSAL=="gene")
#		{
#			glabels = 1:(nGene+1)
#			glabels = glabels[-i]
#			i.new = sample (glabels, size=1)
#			ci.new = ifelse (i.new>nGene, nClust+1, c.curr[i.new])
#		}
#		else
#		{
			clabels = 1:(nClust+1)
			clabels = clabels[-c.curr[i]]
			ci.new = sample (clabels, size=1)
#		}

		if (PRINT) cat (ci.new, "\n")
		
		if (ci.new!=c.curr[i])
		{
			c.new = c.curr
			if (ci.new==nClust+1)
			{
				# cases 1 & 2: proposal=nClust+1, do following:
				# generate new cluster-specific parameter from prior
				# if current cluster is size 1, then partition ratio is alpha
				# if current cluster has size>1, then partition ratio is alpha/(current size-1)
#				logPartRatio = log (ifelse (c.counts[c.curr[i]]>1, alpha.curr/(c.counts[c.curr[i]]-1), 1))
				if (c.counts[c.curr[i]]>1)
				{
					logPartRatio = log (alpha.curr / (c.counts[c.curr[i]]-1))
					logPropRatio = 0
				}
				else
				{
					logPartRatio = 0
					logPropRatio = log (nClust) - log (nClust+1)
				}
				par.i.new = simulateFromPrior (par.prior, times, PRIOR.MODEL)
			}
			else
			{
				# cases 3 & 4: proposal<=nClust and different from current value,
				# find size of parameter of proposal cluster
				# compute partition ratio
#				logPartRatio = log (ifelse (c.counts[c.curr[i]]>1, c.counts[ci.new]/(c.counts[c.curr[i]]-1), c.counts[ci.new]/alpha.curr))
				if (c.counts[c.curr[i]]>1)
				{
					logPartRatio = log (c.counts[ci.new]/(c.counts[c.curr[i]]-1))
					logPropRatio = 0
				}
				else
				{
					logPartRatio = log (c.counts[ci.new]/alpha.curr)
					logPropRatio = log (nClust) - log (nClust-1)
				}
				par.i.new = par.unique.curr[ci.new, ]
			}

			# compute log Hastings ratio and log acceptance probability
			logDataRatio = computeClusterParHastingsRatio (par.i.new, par.unique.curr[c.curr[i],], as.vector (data[i,]), nTime, nRep)
			logHastingsRatio = logDataRatio + logPartRatio + logPropRatio
			accept.prob.log = min (0, logHastingsRatio)
				
			# decide whether to accept c*_i or not
			if (log (runif (1)) < accept.prob.log)
			{
				c.new[i] = ci.new
				if (ci.new == nClust+1)
					par.unique.curr = rbind (par.unique.curr, par.i.new)
			}
				
			if (PRINT) cat ("i=", i, ":", c(c.new, par.i.new), "\n")
				
			c.curr = c.new
		}
		else
			# case 0: proposal same as current value, no action
			need.update = 0
		
	}
		
	if (PRINT) cat ("\n")
	
	return (list (c.curr=c.curr, par.unique.curr=par.unique.curr))
}

