stageTwoAnalysis <-
function(stageOne, explanatory, response, type="IR-wald", level=0.95, combineData=FALSE) {
	cl1<-match.call(expand.dots=TRUE)
	Y_0 <- stageOne$threshold
	C_1 <- stageOne$C_1
	gamma1=1/3
	if (combineData) {
		explanatory <- c(explanatory , stageOne$X1[stageOne$X1 > stageOne$L1 & stageOne$X1 < stageOne$U1])
		response <- c(response , stageOne$Y1[stageOne$X1 > stageOne$L1 & stageOne$X1 < stageOne$U1])		
	}
	if (type == "IR-wald") {
		CI <- waldConfidenceInterval_ir_stageTwo(explanatory, response, Y_0, level=level, gamma1 <- gamma1, C_1 <- C_1, n1=length(stageOne$X1))
		return(structure(list(L2=CI$lower, U2=CI$upper, estimate=CI$estimate, threshold=Y_0, level=level, X1=stageOne$X1, Y1=stageOne$Y1, X2=explanatory, Y2=response, L1=stageOne$L1, U1=stageOne$U1, call=cl1, sigmaSq=CI$sigmaSq, deriv_d0=CI$deriv_d0), class="twostageTE"))
	}
	else if (type == "IR-likelihood") {
		CI <- likelihoodConfidenceInterval(explanatory, response, Y_0, level=level)
		return(structure(list(L2=CI$lower, U2=CI$upper, estimate=CI$estimate, threshold=Y_0, level=level, X1=stageOne$X1, Y1=stageOne$Y1, X2=explanatory, Y2=response, L1=stageOne$L1, U1=stageOne$U1, call=cl1, sigmaSq=CI$sigmaSq, deriv_d0=CI$deriv_d0), class="twostageTE"))
	}
#	else if (type == "SIR") {
#		CI <- waldConfidenceInterval_sir_stageTwo(explanatory=explanatory, response=response, Y_0=Y_0, gamma1 <- gamma1, C_1 <- C_1, level=level)
#		return(structure(list(L2=CI$lower, U2=CI$upper, estimate=CI$estimate, threshold=Y_0, level=level, X1=stageOne$X1, Y1=stageOne$Y1, X2=explanatory, Y2=response, L1=stageOne$L1, U1=stageOne$U1,call=cl1, sigmaSq=CI$sigmaSq, deriv_d0 <- CI$deriv_d0), class="twostageTE"))
#	}
	else if (type == "locLinear") {
		CI <- linearBootstrapConfidenceInterval_stageTwo(explanatory=explanatory, response=response, Y_0=Y_0, level=level)
		return(structure(list(L2=CI$lower, U2=CI$upper, estimate=CI$estimate, threshold=Y_0, level=level, X1=stageOne$X1, Y1=stageOne$Y1, X2=explanatory, Y2=response, L1=stageOne$L1, U1=stageOne$U1,call=cl1, sigmaSq=CI$sigmaSq, deriv_d0=CI$deriv_d0), class="twostageTE"))
	}
#	else if (type == "IR-linear") {
#		CI <- linearConfidenceInterval_stageTwo(explanatory, response, Y_0, level=level)
#		if (verbose) 
#			cat(sprintf("Second stage with local linear, n2=%d samples yields CI: [%.3f,%.3f], with estimate=%.3f\n",length(response),CI$lower,CI$upper, CI$estimate))
#		return(list(L=CI$lower, U=CI$upper, estimate=CI$estimate))
#	}
	else cat("stageOneAnalysis: type should be either 'IR-wald','IR-likelihood' or 'locLinear'")
}
