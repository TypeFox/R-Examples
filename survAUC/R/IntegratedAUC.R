

##############################################################
## Integrated AUC
##############################################################
# AUC	- vector of AUC´s
# times	- vector of times
# S		- vector of survival probability
# tmax	- maximum timepoint



IntAUC <- function( AUC, times, S, tmax, auc.type="cumulative")
{
	n_S <- length(S)
	n_AUC <- length(AUC)
	n_times <- length(times)
	if(!((n_S == n_AUC) && (n_AUC == n_times)))
		stop("AUC, times and S must be the same length!")
	auc.type <- charmatch( auc.type, c("cumulative","incident") )
	if (is.na(auc.type))
		stop("auc.type must be one of 'cumulative' or 'incident'")
	maxI <- sum( times <= tmax )
	ind_S <- S[min(maxI+1,length(S))]
	iAUC <- .C("int_auc",0.0,
			   as.numeric(AUC),
			   as.numeric(times),
			   as.numeric(S),
			   as.numeric(tmax),
			   as.integer(n_S),
			   as.integer(maxI),
			   as.numeric(ind_S),
			   as.integer(auc.type-1),
			   PACKAGE="survAUC")
	iAUC[[1]]
}

