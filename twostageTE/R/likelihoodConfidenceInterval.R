likelihoodConfidenceInterval <-
function(explanatory, response, Y_0, level=NA) {
	if (is.na(level))
		level <- 0.95
	RVforLR_realizations <- NULL; rm(RVforLR_realizations); # Dummy to trick R CMD check 
	data("RVforLR_realizations", envir =environment())
	D <- quantile(RVforLR_realizations, level)
	n <- length(response)
	ind <- order(explanatory, decreasing=FALSE)
	if (sum(diff(ind) < 0) != 0) {
		explanatory <- explanatory[ind]
		response <- response[ind]
	}
	fit <- threshold_estimate_ir(explanatory, response, Y_0)
	sigmaSq <- estimateSigmaSq(explanatory, response)$sigmaSq

	likelihoodRatio <- 
	function(explanatory, response, X_0, Y_0, sigmaSq) {
		logLikelihood <-
		function (Y, Y_hat) {
			-1/(2*sigmaSq)*sum( (Y - Y_hat)^2)
		}

		unconstrainedLikelihood <-
		function(explanatory, response) {  
		#	if (is.na(sigmaSq))
	#			sigmaSq <- estimateSigmaSq(explanatory,response)$sigmaSq
			fit <- pava(explanatory, response)
			tmp <- logLikelihood(fit$response_obs, fit$y)
			return(list(x=fit$x,y_hat=fit$y,y=fit$response_obs, logLikelihood=tmp))
		}
	
		constrainedLikelihood <-
		function(explanatory, response, X_0, Y_0) {
		#	if (is.na(sigmaSq))
	#			sigmaSq <- estimateSigmaSq(explanatory,response)$sigmaSq
			fit <- pava(explanatory, response, X_0, Y_0)
			tmp <- logLikelihood(fit$response_obs, fit$y)
			return(list(x=fit$x,y_hat=fit$y,y=fit$response_obs, logLikelihood=tmp))
		}
	

		unconst <- unconstrainedLikelihood(explanatory,response)
		const <- constrainedLikelihood(explanatory,response, X_0, Y_0)
		return(unconst$logLikelihood - const$logLikelihood)
	}

	
	## These intervals will be slightly conservative, since we are walking along the x-values and stop at the first index such that likelihoodRatio(x,y) >= D
	i =fit$index+1
	lrt_tmp <- 0
	while (lrt_tmp < D && i < n) { 
		lrt_tmp <- likelihoodRatio(explanatory, response, explanatory[i], Y_0, sigmaSq)		
		i <- i + 1
	}
	right <- explanatory[min(i,n)]
	i =fit$index-1
	lrt_tmp <- 0
	while (lrt_tmp < D && i > 1) {
		lrt_tmp <- likelihoodRatio(explanatory, response, explanatory[i], Y_0, sigmaSq)		
		i <- i - 1
	}
	left <- explanatory[max(i,1)]

	return(list(estimate=fit$threshold_estimate_explanatory,lower=left, upper=right, sigmaSq=sigmaSq, deriv_d0=NA ))
}
