stageOneAnalysis <-
function(explanatory, response, threshold, type="IR-wald", level=0.99) {
	cl1<-match.call(expand.dots=TRUE)
	if (type == "IR-wald") {
		CI <- waldConfidenceInterval_ir_stageOne(explanatory, response, threshold, level=level)
		return(structure(list(L1=CI$lower, U1=CI$upper, estimate=CI$estimate, C_1=CI$C_1, threshold=threshold, level=level, X1=explanatory, Y1=response, X2=NA, Y2=NA, L2=NA, U2=NA, call=cl1, sigmaSq=CI$sigmaSq, deriv_d0=CI$deriv_d0), class="twostageTE"))
	}
#	else if (type == "IR-wald-linear") {
#		CI <- waldConfidenceInterval_ir_linear_stageOne(explanatory, response, Y_0, level=level)
#		if (verbose) 
#			cat(sprintf("First stage with IR, n1=%d samples yields CI: [%.3f,%.3f], with estimate=%.3f\n",length(response),CI$lower,CI$upper, CI$estimate))
#		return(list(L=CI$lower, U=CI$upper, estimate=CI$estimate))
#	}
	else if (type == "IR-likelihood") {
		CI <- likelihoodConfidenceInterval(explanatory, response, threshold, level=level)
		return(structure(list(L1=CI$lower, U1=CI$upper, estimate=CI$estimate, threshold=threshold, level=level, X1=explanatory, Y1=response, X2=NA, Y2=NA, L2=NA, U2=NA, call=cl1, sigmaSq=CI$sigmaSq, deriv_d0=CI$deriv_d0), class="twostageTE"))
	}
#else if (type == "SIR") {
#		CI <- waldConfidenceInterval_sir_stageOne(explanatory, response, threshold, level=level)
#		return(structure(list(L1=CI$lower, U1=CI$upper, estimate=CI$estimate, threshold=threshold, level=level, X1=explanatory, Y1=response, X2=NA, Y2=NA, L2=NA, U2=NA, call=cl1, sigmaSq=CI$sigmaSq, deriv_d0 <- CI$deriv_d0), class="twostageTE"))
#	}
	else cat("stageOneAnalysis: type should be either 'IR-wald','IR-likelihood'")
}
