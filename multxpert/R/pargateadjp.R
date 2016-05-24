# ParGateAdjP function computes adjusted p-values and generates decision rules
# for multistage parallel gatekeeping procedures in hypothesis testing problems
# with multiple families of null hypotheses (null hypotheses are assumed
# to be equally weighted within each family)
pargateadjp<-function(gateproc, independence, alpha=0.05, printDecisionRules=FALSE)
# GATEPROC, List of gatekeeping procedure parameters
# INDEPENDENCE, Boolean indicator (TRUE, Independence condition is imposed; FALSE,
# Independence condition is not imposed)
# ALPHA: Global familywise error rate
# PRINTDECISIONRULES: Boolean indicator which controls printing of decision rules
{
	# Number of families
	nfams<-length(gateproc)
	if (nfams<=1) stop("Function requires more than one family of null hypotheses")

	for (i in 1:nfams)
	{
		pr<-gateproc[[i]]$proc
		if (pr!="Bonferroni" & pr!="Holm" & pr!="Hommel" & pr!="Hochberg" & pr!="Fallback")
			stop("Procedure name is not recognized. ParGateAdjP function supports only the Bonferroni, Holm, Hommel, Hochberg and fallback procedures")
	}

	if (alpha <= 0) stop("Alpha must be positive")
	if (alpha >= 1) stop("Alpha must be less than 1")

	temp<-gateproc

	for (i in 1:nfams)
	{
		# Number of null hypotheses
		nhyps<-length(temp[[i]]$rawp)
		adjp<-rep(0,nhyps)
		# Placeholder for adjusted p-values
		gateproc[[i]][5]<-list(adjp=adjp)

		for (j in 1:nhyps)
		{

			# Find the lowest alpha level at which the current null hypothesis is rejected
			upper<-1
			lower<-0
			for (k in 1:20)
			{
				current<-(lower+upper)/2
				# Evaluate decision rules for the multistage parallel gatekeeping procedure
				res<-pargateeval(temp,current,independence)
				# Rejection decision for the current null hypothesis
				if (independence==TRUE | i==nfams) reject<-res[[i]][[7]][j]
				# Rejection decisions after retesting if the independence condition is not imposed
				if (independence==FALSE & i<nfams)
				{
					# If the current null hypothesis was retested
					if (res[[2*nfams-i]][[6]]>0) reject<-res[[2*nfams-i]][[7]][j] else reject<-res[[i]][[7]][j]
				}
				# Update the interval
				if (reject==TRUE) upper<-current
				if (reject==FALSE) lower<-current
			}

			# Global adjusted p-value
			gateproc[[i]][[5]][j]<-(lower+upper)/2
		}
	}

	# Build a data frame with the raw and global adjusted p-values
	count<-0
	for (i in 1:nfams) {
		count<-count+length(gateproc[[i]]$rawp)
	}
        result <- data.frame()
	k<-1
	for (i in 1:nfams)
	{
		# Number of null hypotheses
		nhyps<-length(gateproc[[i]]$rawp)
		for (j in 1:nhyps)
		{
			result[k,1]<-gateproc[[i]]$label
			result[k,2]<-gateproc[[i]]$proc
			result[k,3]<-gateproc[[i]]$procpar
			result[k,4]<-round(gateproc[[i]]$rawp[j], 4)
			result[k,5]<-round(gateproc[[i]][[5]][j], 4)
			k<-k+1
		}
	}
	names(result)[1]<-"Family"
	names(result)[2]<-"Procedure"
	names(result)[3]<-"Parameter"
	names(result)[4]<-"Raw.pvalue"
	names(result)[5]<-"Adj.pvalue"

	if(printDecisionRules==TRUE) { pargaterule(gateproc,alpha,independence)}

	return(result=result)
}
# End of pargateadjp
