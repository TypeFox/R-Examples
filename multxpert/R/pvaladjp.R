# PValAdjP function computes adjusted p-values and generates decision rules
# for the Bonferroni, Holm, Hommel, Hochberg, fixed-sequence and fallback procedures
pvaladjp<-function(rawp,weight=rep(1/length(rawp), length(rawp)),alpha=0.05,
		           proc=c("Bonferroni", "Holm", "Hommel", "Hochberg", "Fixed-sequence", "Fallback"),
				   printDecisionRules=FALSE)
# RAWP, Vector of raw p-values
# WEIGHT, Vector of hypothesis weights
# ALPHA, Familywise error rate
# PROC, Procedure name
# PRINTDECISIONRULES: Boolean indicator which controls printing of decision rules

{

	# Number of null hypotheses
	m<-length(rawp)

	if (m==0) stop("No p-values are specified")

	for (i in 1:m)
	{
		if (rawp[i]<0) stop("P-values must be positive")
		if (rawp[i]>1) stop("P-values must be less than 1")
	}

	index <-order(rawp)

	if (alpha <= 0) stop("Alpha must be positive")
	if (alpha >= 1) stop("Alpha must be less than 1")

	# Number of weights
	nweis<-length(weight)

	if (m!=nweis) stop("RAWP and WEIGHT vectors have different lengths")

	if (sum(weight)>1) stop("Sum of hypothesis weights must be <=1")

	for (i in 1:nweis)
	{
		if (weight[i]<0) stop("Hypothesis weights must be >=0")
	}

	if(!all(proc %in% c("Bonferroni", "Holm", "Hommel", "Hochberg", "Fixed-sequence", "Fallback")))
		stop("Procedure name is not recognized. PvalAdjp function supports the Bonferroni, Holm, Hommel, Hochberg, Fixed-sequence, and Fallback procedures")

	# number of procedures specified
	nproc <- length(proc)

	# set up matrix to contain adjusted p-values
	adjp <- matrix(0,m,nproc)
	dimnames(adjp) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

	if (is.element("Bonferroni", proc)) {
		adjp[,"Bonferroni.adj.pvalue"]<-pvaltrunc(rawp,weight,"Holm",0)
	}
	if (is.element("Holm", proc)) {
		adjp[,"Holm.adj.pvalue"]<-pvaltrunc(rawp,weight,"Holm",1)
	}
	if (is.element("Hommel", proc)) {
		adjp[,"Hommel.adj.pvalue"]<-pvaltrunc(rawp,weight,"Hommel",1)
	}
	if (is.element("Hochberg", proc)) {
		adjp[,"Hochberg.adj.pvalue"]<-pvaltrunc(rawp,weight,"Hochberg",1)
	}
	if (is.element("Fixed-sequence", proc)) {
		adjp[,"Fixed-sequence.adj.pvalue"]<-pvaltrunc(rawp,weight,"Fixed-sequence",1)
	}
	if (is.element("Fallback", proc)) {
		adjp[,"Fallback.adj.pvalue"]<-pvaltrunc(rawp,weight,"Fallback",1)
	}

	# Data frame returned by the function
	result<-data.frame(round(rawp,4),round(weight,4), round(adjp,4))
	names(result)[1]<-"Raw.pvalue"
	names(result)[2]<-"Weight"

	if(printDecisionRules==TRUE) {
            if(length(unique(round(weight,3))) > 1) stop("Weights must be equal for decision rule calculations to be valid")
            pvalrule(rawp, alpha, proc)
        }

	return(result=result)
}
# End of pvaladjp
