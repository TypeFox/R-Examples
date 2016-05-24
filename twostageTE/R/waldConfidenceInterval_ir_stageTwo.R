waldConfidenceInterval_ir_stageTwo <-
function(explanatory, response, Y_0, gamma1, C_1, n1, level=NA) {
	if (is.na(level)) {
		level <- 0.95
	}
	alpha <- 1 - level
	## Import previously computed Chernoff quantiles, provided by Groeneboom and Wellner
#	realizations <- read.table("results_chernoff", header=TRUE)
	chernoff_realizations <- NULL; rm(chernoff_realizations); # Dummy to trick R CMD check 
	data("chernoff_realizations", envir =environment())

	ind <- min(which(chernoff_realizations$DF - (1-alpha/2) >= 0))
	q <- chernoff_realizations$xcoor[ind]
	n <- length(response)
	
	fit <- threshold_estimate_ir(explanatory, response, Y_0)

	phi_0 <- C_1*n1*(n^(-1)) # densTMP # Phi(0) on page 9

	sigmaSq <- estimateSigmaSq(explanatory, response)$sigmaSq
	deriv_d0 <- estimateDeriv(explanatory, response, fit$threshold_estimate_explanatory, sigmaSq) 
	C_di <- (4*sigmaSq / (deriv_d0^2) )^(1.0/3.0)

	n <- length(explanatory)
	p <- gamma1/(1.0+gamma1)
#	C_di <- C_1 / q_s1
	C_di2 <- C_di * (C_1 / ((1-p)*p^(gamma1) * phi_0))#^(1/3)
	
	band <- n^(-1*(1 + gamma1)/3.0) * C_di2 * q
	return(list(estimate=fit$threshold_estimate_explanatory,lower=max(min(explanatory),fit$threshold_estimate_explanatory - band), upper=min(max(explanatory),fit$threshold_estimate_explanatory + band), sigmaSq=sigmaSq, deriv_d0=deriv_d0))
}
