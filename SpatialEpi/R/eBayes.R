eBayes <-
function(Y, E, Xmat=NULL){
	# Check for covariates
	if (is.null(Xmat)) {
		mod <- glm.nb(Y ~ 1+offset(log(E)), link=log)
	}else{
		mod <- glm.nb(Y ~ Xmat+offset(log(E)), link=log)
	}

	# Get MLE's of parameters
	alpha <- mod$theta
	muhat <- mod$fitted/E
	
	# Generate weights between global and local SMR
	wgt <- E*muhat/(alpha+E*muhat); SMR <- Y/E
	RR <- as.numeric(wgt*SMR + (1-wgt)*muhat)
	RRmed <- qgamma(0.5,alpha+Y,(alpha+E*muhat)/muhat)
	
	# Output results
	list(RR=RR, RRmed=RRmed, beta=mod$coeff, alpha=alpha, SMR=SMR)
}
