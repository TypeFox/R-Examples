threshold_estimate_ir <-
function (explanatory, response, Y_0) {
	n <- length(response)
	if (sum(response < Y_0) == n ) {
		warning("Y_0 is outside (above) observed data")
		list(threshold_estimate_explanatory=max(explanatory), threshold_estimate_response=max(response), threshold=Y_0, Y_hat=max(response), index=n)
	}
	else if (sum(response >= Y_0) == n ) {
		warning("Y_0 is outside (below) observed data")
		list(threshold_estimate_explanatory=min(explanatory), threshold_estimate_response=min(response), threshold=Y_0, Y_hat=min(response), index=1)
	}	
	else {
		fit <- pava(explanatory, response)
		if (sum(fit$y >= Y_0) == 0) {
			warning("Threshold estimate is on the lower boundary point of the estimated isotonic regression function")
			ind <- n
			estim_x <- fit$x[ind]
		}
		else if (sum(fit$y <= Y_0) == 0) {
			warning("Threshold estimate is on the upper boundary point of the estimated isotonic regression function")
			ind <- min(which(fit$y >= Y_0))  ## fit$y is already sorted
			estim_x <- fit$x[ind]
		}
		else {
			ind <- min(which(fit$y >= Y_0))  ## fit$y is already sorted
			estim_x <- fit$x[ind]
		}
#cat(sprintf("\nd_0=%.2f, m(d_0)=%.2f, est(d_0)=%.2f", sqrt(Y_0), Y_0, estim_x))
		list(threshold_estimate_explanatory = estim_x, threshold_estimate_response = fit$y[ind], threshold = Y_0, Y_hat = fit$y, index=ind)
	}
}
