# Bonf function is a secondary function which computes the weighted Bonferroni
# p-value for an intersection hypothesis
bonf<-function(p,w)
# P, Vector of raw p-values
# W, Vector of hypothesis weights
{
	# Number of null hypotheses included in the intersection
	k<-length(w[w!=0])
	if (k>0) bonf<-min(p[w!=0]/w[w!=0])
	if (k==0) bonf<-1
	return(bonf)
}
# End of bonf




# Simes function is a secondary function which computes a truncated version of the
# weighted Simes p-value for an intersection hypothesis
simes<-function(p,w,gamma)
# P, Vector of raw p-values
# W, Vector of hypothesis weights
# GAMMA, Truncation parameter
{
	# Number of null hypotheses included in the intersection
	k<-length(w[w!=0])
	if (k>1)
	{
		temp<-matrix(0,3,k)# TODO: Add comment
		# P-values
		temp[1,]<-p[w!=0]
		# Weights
		temp[2,]<-w[w!=0]
		# Normalized weights
		temp[3,]<-w[w!=0]/sum(w[w!=0])
		# Sort by p-values
		sorted<-temp[,order(temp[1,])]
		numer<-sorted[1,]
		denom<-gamma*cumsum(sorted[3,])+(1-gamma)*sorted[2,]
		simes<-min(numer/denom)
	}
	if (k==1)
	{
		numer<-p[w!=0]
		denom<-gamma+(1-gamma)*w[w!=0]
		simes<-numer/denom
	}
	if (k==0) simes<-1
	return(simes)
}
# End of simes

# IncSimes function is a secondary function which computes a truncated version of the weighted
# incomplete Simes p-value for an intersection hypothesis
incsimes<-function(p,w,gamma)
# P, Vector of raw p-values
# W, Vector of hypothesis weights
# GAMMA, Truncation parameter
{
	# Number of null hypotheses included in the intersection
	k<-length(w[w!=0])
	if (k>1)
	{
		temp<-matrix(0,3,k)
		# P-values
		temp[1,]<-p[w!=0]
		# Weights
		temp[2,]<-w[w!=0]
		# Normalized weights
		temp[3,]<-w[w!=0]/sum(w[w!=0])
		# Sort by p-values
		sorted<-temp[,order(temp[1,])]
		modw<-w[w!=0]
		modw[1]<-0
		modw[2:k]<-sorted[3,1:k-1]
		numer<-sorted[1,]
		denom<-gamma*sorted[3,]/(1-cumsum(modw))+(1-gamma)*sorted[2,]
		incsimes<-min(numer/denom)
	}
	if (k==1)
	{
		numer<-p[w!=0]
		denom<-gamma+(1-gamma)*w[w!=0]
		incsimes<-numer/denom
	}
	if (k==0) incsimes<-1
	return(incsimes)
}
# End of incsimes

# PValTrunc function is a secondary function which computes adjusted p-values for
# truncated versions of the Holm, Hommel, Hochberg, fixed-sequence and fallback procedures
# in general hypothesis testing problems (equally or unequally weighted null hypotheses)
pvaltrunc<-function(rawp,weight,proc,gamma)
# RAWP, Vector of raw p-values
# WEIGHT, Vector of hypothesis weights
# PROC, Procedure name
# GAMMA, Truncation parameter
{
	# Number of null hypotheses
	nhyps<-length(rawp)

	# Number of weights
	nweis<-length(weight)

	# Number of intersection hypotheses
	nints<-2**nhyps-1

	# Decision matrix
	h<-matrix(0,nints,nhyps)
	for(i in 1:nhyps)
	{
		for(j in 0:(nints-1))
		{
			k<-floor(j/2**(nhyps-i))
			if (k/2==floor(k/2)) h[j+1,i]<-1
		}
	}

	# Temporary vector of weights
	tempw<-matrix(0,1,nhyps)

	# Holm procedure
	if (proc=="Holm")
	{
		# Local p-values for Holm procedure
		holmp<-matrix(0,nints,nhyps)

		for(i in 1:nints)
		{
			tempw<-(weight*h[i,]/sum(weight*h[i,]))*gamma+(1-gamma)*weight*h[i,]
			holmp[i,]<-h[i,]*bonf(rawp,tempw)
		}

		# Compute adjusted p-values
		holm<-rep(0,nhyps)
		for(i in 1:nhyps) holm[i]<-pmin(1, max(holmp[,i]))

		# Delete temporary objects
		rm(holmp)

		# List of adjusted p-values
		res<-holm
	}

	# Hommel procedure
	else if (proc=="Hommel")
	{
		# Local p-values for Hommel procedure
		hommp<-matrix(0,nints,nhyps)

		for(i in 1:nints) hommp[i,]<-h[i,]*simes(rawp,weight*h[i,],gamma)

		# Compute adjusted p-values
		hommel<-rep(0,nhyps)
		for(i in 1:nhyps) hommel[i]<-max(hommp[,i])

		# Delete temporary objects
		rm(hommp)

		# List of adjusted p-values
		res<-hommel
	}

	# Hochberg procedure
	else if (proc=="Hochberg")
	{
		# Local p-values for Hochberg procedure
		hochp<-matrix(0,nints,nhyps)

		for(i in 1:nints) hochp[i,]<-h[i,]*incsimes(rawp,weight*h[i,],gamma)

		# Compute adjusted p-values
		hochberg<-rep(0,nhyps)
		for(i in 1:nhyps) hochberg[i]<-max(hochp[,i])

		# Delete temporary objects
		rm(hochp)

		# List of adjusted p-values
		res<-hochberg
	}

	# Fallback procedure
	else if (proc=="Fallback")
	{
		# Local p-values for fallback procedure
		fallp<-matrix(0,nints,nhyps)
		fallwin<-matrix(0,nhyps,nhyps)

		# Compute window variables
		for(i in 1:nhyps)
		{
			for(j in 1:nhyps)
			{
				if (i>=j) fallwin[i,j]<-1
			}
		}

		for(i in 1:nints)
		{
			tempw[]<-0
			for(j in 1:nhyps)
			{
				if (h[i,j]==1) tempw[j]<-(sum(weight*fallwin[j,])-sum(tempw))
				if (h[i,j]==0) tempw[j]<-0
			}
			tempw<-tempw*gamma+(1-gamma)*weight*h[i,]
			fallp[i,]<-h[i,]*bonf(rawp,tempw)
		}

		# Compute adjusted p-values
		fallback<-rep(0,nhyps)
		for(i in 1:nhyps) fallback[i]<-pmin(1,max(fallp[,i]))

		# Delete temporary objects
		rm(fallp)
		rm(fallwin)

		# List of adjusted p-values
		res<-fallback
	}

	# Fixed-sequence procedure
	else if (proc=="Fixed-sequence")
	{
		# Compute adjusted p-values using a recursive algorithm
		fixedseq<-rep(0,nhyps)
		fixedseq[1]<-rawp[1]
		if (nhyps>1)
		{
			for(i in 2:nhyps) fixedseq[i]<-max(fixedseq[i-1],rawp[i])
		}

		# List of adjusted p-values
		res<-fixedseq
	}

	rm(tempw)
	return(res)
}
# End of pvaltrunc
