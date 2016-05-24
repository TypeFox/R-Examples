untrans.theta <-
function(theta, scale)
{

	ind <- which(scale == "log")
	if (length(ind)) 
	theta[ind] <- exp(theta[ind])
	ind <- which(scale == "logistic")
	if (length(ind)) 
	theta[ind] <- exp(theta[ind])/(1 + exp(theta[ind]))
	ind <- which(scale == "logistic180")
	if (length(ind)) 
		theta[ind] <- 180 * exp(theta[ind])/(1 + exp(theta[ind]))
	theta

}

