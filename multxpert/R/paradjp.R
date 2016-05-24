# ParAdjP function computes adjusted p-values for the single-step Dunnett
# procedure and step-down Dunnett procedure in one-sided hypothesis testing
# problems with a balanced one-way layout and equally weighted null hypotheses
paradjp<-function(stat,n,proc=c("Single-step Dunnett", "Step-down Dunnett"))
# STAT, Vector of test statistics
# N, Common sample size in each treatment group
# PROC, Procedure name
{

	if (n<=0) stop("Sample size must be positive")

	# Number of null hypotheses
	m<-length(stat)

        if(m==0) stop("No test statistics are specified")

        if(length(n)==0) stop("No sample size specified")

	# Number of degrees of freedom
	nu<-(m+1)*(n-1)
	# Raw p-values
	rawp<-1-pt(stat,2*(n-1))
	# Adjusted p-values
	adjp<-rep(0,m)

        if(!all(proc %in% c("Single-step Dunnett", "Step-down Dunnett"))) stop("Procedure name is not recognized. ParAdjP function supports only the single-step Dunnett and step-down Dunnett procedures")

        # Number of procedures specified
        nproc <- length(proc)

        # Set up matrix to contain adjusted p-values
        adjp <- matrix(0,m,nproc)
        dimnames(adjp) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

	if (is.element("Single-step Dunnett", proc))
        {
        #for (i in 1:m) adjp[i]<-1-pdunnett(stat[i],nu,m)
        adjp[,"Single-step Dunnett.adj.pvalue"] <- sapply(stat, function(x) {1-pdunnett(x,nu,m)})
        }
	if (is.element("Step-down Dunnett", proc))
        {
        adjptmp <- rep(NA, length(stat))
		# Sort test statistics from largest to smallest
		or<-order(stat,decreasing=TRUE)
		sorted<-stat[or]
		for (i in 1:m)
		{
			if (i==1)
			{
				#adjp[1,"Step-down Dunnett.adj.pvalue"] <- 1-pdunnett(sorted[1],nu,m)
				adjptmp[1] <- 1-pdunnett(sorted[1],nu,m)
                                #print(1-pdunnett(sorted[1],nu,m))
				#maxp<-adjp[1, "Step-down Dunnett.adj.pvalue"]
				maxp<-adjptmp[1]
			}
			if (i>1 & i<m)
			{
				#adjp[i, "Step-down Dunnett.adj.pvalue"] <- max(maxp,1-pdunnett(sorted[i],nu,m-i+1))
				adjptmp[i] <- max(maxp,1-pdunnett(sorted[i],nu,m-i+1))

				#maxp<-max(maxp,adjp[i, "Step-down Dunnett.adj.pvalue"])
				maxp<-max(maxp,adjptmp[i])
			}
			if (i==m) {
                            #adjp[m, "Step-down Dunnett.adj.pvalue"] <- max(maxp,1-pt(sorted[m],nu))
                            adjptmp[m] <- max(maxp,1-pt(sorted[m],nu))
                        }
		}
		# Return to original ordering
		temp<-adjptmp
		adjptmp[or]<-temp
                adjp[,"Step-down Dunnett.adj.pvalue"] <- adjptmp
	}

	# Data frame returned by the function
	result<-data.frame(stat, round(rawp, 4), adjp)
	names(result)[1]<-"Test.statistic"
	names(result)[2]<-"Raw.pvalue"
	#names(result)[3]<-"Adj.pvalue"

	return(result=result)
}
# End of paradjp
