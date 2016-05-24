# The PvalCI function computes one-sided multiplicity-adjusted confidence
# intervals (simultaneous confidence intervals) for the Bonferroni, Holm and fixed-sequence
# procedures in general hypothesis testing problems with equally or
# unequally weighted null hypotheses
pvalci<-function(rawp,est,stderror,weight=rep(1/length(rawp),length(rawp)),covprob=0.975,proc=c("Bonferroni", "Holm", "Fixed-sequence"))
# RAWP, Vector of raw p-values
# EST, Vector of point estimates
# STDERROR, Vector of standard errors associated with the point estimates
# WEIGHT, Vector of hypothesis weights
# COVPROB, Simultaneous coverage probability
# PROC, Procedure name
    {

    # Number of null hypotheses
    m<-length(rawp)

    if (m==0) stop("No p-values are specified")

    for (i in 1:m)
        {
        if (rawp[i]<0) stop("P-values must be positive")
        if (rawp[i]>1) stop("P-values must be less than 1")
        }

    if (m!=length(weight)) stop("RAWP and WEIGHT vectors have different lengths")
    if (m!=length(est)) stop("RAWP and EST vectors have different lengths")
    if (m!=length(stderror)) stop("RAWP and STDERROR vectors have different lengths")

    if (sum(weight)>1) stop("Sum of hypothesis weights must be <=1")

    for (i in 1:length(weight))
        {
        if (weight[i]<0) stop("Hypothesis weights must be >=0")
        }

    if (covprob>=1) stop("Simultaneous coverage probability must be <1")
    if (covprob<=0) stop("Simultaneous coverage probability must be >0")

    if (!all(proc %in% c("Bonferroni", "Holm", "Fixed-sequence")))
        stop("Procedure name is not recognized")
    #if (proc!="Bonferroni" & proc!="Holm" & proc!="Fixed-sequence") stop("Procedure name is not recognized")

    # number of procedures specified
    nproc <- length(proc)

        # set up matrix to contain confidence limits
    cimat <- matrix(0,m,nproc)
    dimnames(cimat) <- list(NULL, paste(proc, ".conf.limit", sep=""))

 	# One-sided familywise error rate
 	alpha<-1-covprob

    # Compute adjusted p-values
    result <- pvaladjp(rawp=rawp,weight=weight,alpha=alpha,proc=proc)
    adjpmat <- result[, grep(".adj.pvalue", names(result), value=TRUE)]
    #adjp<-out[c(-1,-2)]
    #print(out)

        # Rejection/acceptance of null hypotheses
	#reject<-(adjp<=alpha)

    # Vectors of confidence limits
    ci<-rep(0,m)

 	zero<-rep(0,m)


   # adjp<-out[,3]

    # Bonferroni procedure
	if (is.element("Bonferroni", proc)) {
           reject <- (result[,"Bonferroni.adj.pvalue"] <= alpha)
           cimat[,"Bonferroni.conf.limit"] <-est-(stderror*qnorm(1-(alpha*weight)))
        }

    # Holm procedure
	if (is.element("Holm", proc)) {
            reject <- (result[,"Holm.adj.pvalue"] <= alpha)
        bonfci<-est-(stderror*qnorm(1-(alpha*weight)))
        # All null hypotheses are rejected
 		if (sum(reject)==m) cimat[,"Holm.conf.limit"] <-pmax(zero,bonfci)
        # Some null hypotheses are accepted
        if (sum(reject)<m)
            {
            for(i in 1:m)
                {
                if (reject[i]==1) cimat[i, "Holm.conf.limit"]<-0
                if (reject[i]==0)
                    {
                    adjalpha<-(alpha*weight[i])/sum(weight[reject==0])
                    cimat[i,"Holm.conf.limit"]<-est[i]-(stderror[i]*qnorm(1-adjalpha))
                    }
                }
            }
        }

   	# Fixed-sequence procedure
	if (is.element("Fixed-sequence", proc)) {
            reject <- (result[,"Fixed.sequence.adj.pvalue"] <= alpha)
		# All null hypotheses are accepted
  		if (sum(reject)==0)
            {
            cimat[1,"Fixed-sequence.conf.limit"] <-est[1]-stderror[1]*qnorm(1-alpha)
            for(i in 2:m) cimat[i, "Fixed-sequence.conf.limit"]<-NA
            }
        # All null hypotheses are rejected
		if (sum(reject)==m)
            {
            temp1<-est-stderror*qnorm(1-alpha)
            cimat[,"Fixed-sequence.conf.limit"] <-min(temp1)
            }
        # Some null hypotheses are accepted and some are rejected
        if (sum(reject)>0 & sum(reject)<m)
            {
            cimat[1, "Fixed-sequence.conf.limit"]<-0
            for(i in 2:m)
                {
                if (reject[i]==1) cimat[i, "Fixed-sequence.conf.limit"]<-0
                if (reject[i]==0 & reject[i-1]==1) cimat[i, "Fixed-sequence.conf.limit"]<-est[i]-stderror[i]*qnorm(1-alpha)
                if (reject[i]==0 & reject[i-1]==0) cimat[i, "Fixed-sequence.conf.limit"]<-NA
                }
            }
        }

	# Data frame returned by the function
    #result<-data.frame(rawp,est,stderror,weight,adjp,ci)
    result<-data.frame(rawp,est,stderror,weight,adjpmat,cimat)
    names(result)[1]<-"Raw.pvalue"
    names(result)[2]<-"Estimate"
    names(result)[3]<-"Std.error"
    names(result)[4]<-"Weight"
    #names(result)[5]<-"Adj.pvalue"
    #names(result)[6]<-"Conf.limit"

    return(result=result)
}
# End of pvalci
