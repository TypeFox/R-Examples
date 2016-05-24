# 
# Author: JonathanRosenblatt
###############################################################################



pval2qval<- function(pValues, cutoff){
	requireLibrary('fdrtool')
	fdrtool <- get("fdrtool", envir=asNamespace("fdrtool"))
	qvals<-fdrtool(
			pValues,
			statistic= 'pvalue', 
			plot=FALSE,verbose=FALSE)$qval
	
	if (missing(cutoff)) {
		return(list(qValues=qvals))
	}
	return(list(qValues=qvals, rejected= qvals<=cutoff ))		
	
}



pval2locfdr<- function(pValues, cutoff){
	requireLibrary('fdrtool')
	fdrtool <- get("fdrtool", envir=asNamespace("fdrtool"))
	locfdr<-fdrtool(
			pValues,
			statistic= 'pvalue', 
			plot=FALSE,verbose=FALSE)$lfdr
	
	if (missing(cutoff)) {
		return(list(locFDR=locfdr))
	}
	return(list(locFDR=locfdr, rejected= locfdr<=cutoff ))		
}

mutoss.locfdr <- function() { 
	return(new(Class="MutossMethod",
					label="Local FDR (fdr)",
					callFunction="pval2locfdr",
					output=c("locFDR", "rejected"),  
					info=
							"<h2> Name: </h2> Local fdr.\n
							<h3> Also known as: </h3> fdr, empirical posterior probability of the null. \n
							<h3> Error Type: </h3> Motivated by Bayesian considerations. Does not guarantee control of frequentist error types like FWER or FDR.\n
							<h3> Recommended Usage: </h3> Typically used when a massive amount of hypotheses is being tested as in microarray analyses.\n
							<h3> Related procedures: </h3> See FDR methods for similar procedures for frequentist error control.\n
							<h3> References: </h3> \n
							<ul>
							<li> Efron B., Tibshirani R., Storey J. D. and Tusher, V. (2001).<i> Empirical Bayes Analysis of a Microarray Experiment. </i>\n
								 Journal of the American Statistical Association 96(456):1151-1160. </li>
							</ul>",				
					parameters=list(
							pValues=list(type="numeric"),
							cutoff=list(type="numeric", label="Local fdr cutoff for rejection", optional=TRUE))))
}



mutoss.qvalues <- function() { 
	return(new(Class="MutossMethod",
					label="q Values (Fdr)",
					callFunction="pval2qval",
					output=c("qValues", "rejected"),
					info=
							"<h2> Name: </h2> q-Values.\n
							<h3> Also known as: </h3> \n
							<ul> 
								<li> Estimated pFDR</li>\n 
								<li> Estimated Positive FDR </li>\n
								<li> Empirical tail-area posterior probability of the null</li> \n
							<h3> Error Type: </h3> Motivated by Bayesian considerations. Guarantees FDR control only when masses of hypotheses are being tested.\n
							<h3> Recommended Usage: </h3> Typically used when a massive amount of hypotheses is being tested as in microarray analyses.\n
							<h3> Related procedures: </h3> See FDR methods for similar procedures with frequentist error control.\n
							<h3> References: </h3> \n
							<ul>
							<li> Storey, J. D. (2003)<i>The Positive False Discovery Rate: A Bayesian Interpretation and the q-Value.</i>
							 The Annals of Statistics 31(6): 2013-2035. </li>
							</ul>",
					parameters=list(
							pValues=list(type="numeric"),
							cutoff=list(type="numeric", label="q-value (pFDR) cutoff for rejection", optional=TRUE))))
}
