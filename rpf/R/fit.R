# copied from OpenMx R/MxSummary.R

pChiSqFun <- function(x, val, degf, goal){
	goal - pchisq(val, degf, ncp=x)
}

rmseaConfidenceIntervalHelper <- function(chi.squared, df, N, lower, upper){
	# Lower confidence interval
	if( pchisq(chi.squared, df=df, ncp=0) >= upper){ #sic
		lower.lam <- uniroot(f=pChiSqFun, interval=c(1e-10, 1e4), val=chi.squared, degf=df, goal=upper)$root
		# solve pchisq(ch, df=df, ncp=x) == upper for x
	} else{
		lower.lam <- 0
	}
	# Upper confidence interval
	if( pchisq(chi.squared, df=df, ncp=0) >= lower){ #sic
		upper.lam <- uniroot(f=pChiSqFun, interval=c(1e-10, 1e4), val=chi.squared, degf=df, goal=lower)$root
		# solve pchisq(ch, df=df, ncp=x) == lower for x
	} else{
		upper.lam <- 0
	}
	lower.rmsea <- sqrt(lower.lam/(N*df))
	upper.rmsea <- sqrt(upper.lam/(N*df))
	return(c(lower=lower.rmsea, upper=upper.rmsea))
}

computeFitStatistics <- function(likelihood, DoF, numObs,
				 independence, indDoF, saturated=0, satDoF=0) {
	chi <- likelihood - saturated
	chiDoF <- DoF - satDoF # DoF = obsStat-model.ep; satDoF = obsStat-sat.ep; So sat.ep-model.ep == DoF-satDoF
	CFI <- (independence - indDoF - likelihood + DoF)/(independence - indDoF - saturated + satDoF)
	TLI <- 1
	rmseaSquared <- 0
	RMSEA <- 0
	RMSEACI <- c(lower=NA, upper=NA)
	if (!is.na(chiDoF) && chiDoF > 0) {
		TLI <- ((independence-saturated)/(indDoF-satDoF) - (chi)/(DoF-satDoF))/((independence-saturated)/(indDoF-satDoF) - 1)
					# Here we use N in the denominator as given in the original
					# RMSEA paper. The difference between N and N-1 is negligible
					# for sample sizes over 30. RMSEA should not be taken seriously
					# such small samples anyway.
		rmseaSquared <- (chi / (chiDoF) - 1) / numObs
		RMSEACI <- c(lower=NA, upper=NA)
		if (length(rmseaSquared) == 0 || is.na(rmseaSquared) || 
		    is.nan(rmseaSquared)) { 
					# || (rmseaSquared < 0)) { # changed so 'rmseaSquared < 0' yields zero with comment
			RMSEA <- NA
		} else if (rmseaSquared < 0) {
			RMSEA <- 0.0
		} else {
			RMSEA <- sqrt(rmseaSquared)
			ci <- try(rmseaConfidenceIntervalHelper(chi, chiDoF, numObs, .025, .975))
			if (!inherits(ci, "try-error")) RMSEACI <- ci
		}
	}
	list(CFI=CFI, TLI=TLI, RMSEA=RMSEA, RMSEASquared=rmseaSquared, RMSEACI=RMSEACI)
}

catFitStatistics <- function(x) {
	cat("CFI:", x$CFI, '\n')
	cat("TLI:", x$TLI, '\n')
	if (length(x$RMSEASquared) == 1 && !is.na(x$RMSEASquared) && x$RMSEASquared < 0.0) {
		cat("RMSEA: ", x$RMSEA, '*(Non-centrality parameter is negative)', '\n')
	} else {
		cat("RMSEA:  ", x$RMSEA, "  [95% CI (", x$RMSEACI[1], ", ", x$RMSEACI[2], ")]", '\n', sep="")
	}
}
