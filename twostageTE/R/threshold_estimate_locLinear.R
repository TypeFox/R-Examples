threshold_estimate_locLinear <-
function (explanatory, response, Y_0) {
	n <- length(response)
	if (sum(response < Y_0) == n ) {
		list(threshold_estimate_explanatory=max(explanatory), threshold_estimate_response=max(response), threshold=Y_0, Y_hat=max(response), index=n)
	}
	else if (sum(response >= Y_0) == n ) {
		list(threshold_estimate_explanatory=min(explanatory), threshold_estimate_response=min(response), threshold=Y_0, Y_hat=min(response), index=1)
	}	
	else {
		beta <- lm(response ~ explanatory)$coef
		estim_x <- as.numeric((Y_0 - beta[1]) / beta[2])
		list(threshold_estimate_explanatory=estim_x, threshold=Y_0)
	}
}
