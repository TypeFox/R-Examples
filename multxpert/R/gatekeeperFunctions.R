# LevelProc is a secondary function which computes the significance level for 
# the next family for multistage parallel gatekeeping procedures in hypothesis 
# testing problems with equally weighted null hypotheses
levelproc<-function(reject,proc,gamma,level)
# REJECT, Vector of rejection decisions (1, if rejected; 0, if accepted)
# PROC, Procedure name
# GAMMA, Truncation parameter
# LEVEL, Significance level in the current family
{
	# Total number of null hypotheses in the current family
	m<-length(reject)
	# Number of null hypotheses accepted in the current family
	a<-m-sum(reject)
	# Ratio of accepted hypotheses in the current family
	ratio<-a/m
	
	if (ratio==0) error<-0
	if (ratio>0)
	{
		if (proc=="Holm" | proc=="Hommel" | proc=="Hochberg")
		{
			error<-gamma+(1-gamma)*ratio
		}
		else if (proc=="Fallback")
		{
			accept<-1-reject
			sum<-0
			smallest<-1
			for (i in 1:length(accept))
			{
				# Loop over accepted null hypotheses
				if (accept[i]==1)
				{
					# i is not the smallest index
					if (smallest==0)
					{
						largest<-1
						for (j in 1:(i-1))
						{
							if (accept[j]==1) largest<-j
						}
						l<-largest
					}
					
					# i is the smallest index
					if (smallest==1)
					{
						l<-0
						smallest<-0
					}
					sum<-sum+(i-l)
				}
				error<-gamma*(sum/m)+(1-gamma)*ratio
			}
		}
	}
	
	# Significance level in the next family
	nextlevel<-level-error*level
	return(nextlevel)
}
# End of levelproc

# ParGateEval is a secondary function which evaluates decision rules for 
# multistage parallel gatekeeping procedures in hypothesis testing problems 
# with equally weighted null hypotheses
pargateeval<-function(gateproc,alpha,independence)
# GATEPROC, List of gatekeeping procedure parameters
# ALPHA, Global familywise error rate
# INDEPENDENCE, Boolean indicator (TRUE, Independence condition is imposed; FALSE, 
# Independence condition is not imposed)
{
	
	# Number of families
	nfams<-length(gateproc)
	
	# Significance level in the first family
	level<-alpha
	
	# Test families from first to last
	for(i in 1:nfams)
	{
		# Raw p-values in the current family
		p<-gateproc[[i]]$rawp
		# Null hypotheses are equally weighted in the current family
		w<-rep(1/length(p),length(p))
		# Procedure and truncation parameter in the current family
		if (gateproc[[i]]$proc=="Bonferroni") 
		{
			proc<-"Holm"
			gamma<-0
			gateproc[[i]]$proc<-"Holm"
			gateproc[[i]]$procpar<-0
		}
		else
		{
			proc<-gateproc[[i]]$proc
			gamma<-gateproc[[i]]$procpar
		}
		# Compute the local adjusted p-values in the current family
		adjp<-pvaltrunc(p,w,proc,gamma)
		# Rejection decisions in the current family (0 if rejected and 1 if accepted)
		rejection<-(adjp<=level)
		# Save the local adjusted p-values in the current family
		gateproc[[i]][5]<-list(adjp=adjp)
		# Save the significance level in the current family
		gateproc[[i]][6]<-list(level=level)
		# Save the rejection decisions in the current family
		gateproc[[i]][7]<-list(rejection=rejection)
		# Compute the significance level in the next family
		level<-levelproc(rejection,proc,gamma,level)
	}
	
	# Retest families from last to first if the independence condition is not imposed    
	if (independence==FALSE)    
	{        
		for(i in (nfams+1):(2*nfams-1))
		{
			# Family index
			k<-2*nfams-i     
			prevrej<-gateproc[[i-1]][[7]]   
			# Full alpha level if all null hypotheses are rejected in the previous family            
			if (sum(prevrej)==length(prevrej)) level<-alpha else level<-0            
			# Label in the current family
			label<-gateproc[[k]]$label
			# Raw p-values in the current family
			p<-gateproc[[k]]$rawp
			# Null hypotheses are equally weighted in the current family
			w<-rep(1/length(p),length(p))
			# Procedure in the current family
			if (gateproc[[k]]$proc=="Bonferroni") proc<-"Holm" else proc<-gateproc[[k]]$proc
			# Truncation parameter in the current family (regular procedure is used)
			gamma<-1
			# Compute the local adjusted p-values in the current family
			adjp<-pvaltrunc(p,w,proc,gamma)
			# Rejection decisions in the current family (0 if rejected and 1 if accepted)
			rejection<-(adjp<=level)
			# Create a new list for the current family
			gateproc[[i]]<-list(label=label, rawp=p, proc=proc, procpar=gamma, adjp=adjp, level=level, rejection=rejection)            
		}
	}
	
	return(gateproc)
	
}
# End of pargateeval