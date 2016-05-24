estimateSigmaSq <-
function (explanatory, response) {
	tb = table(explanatory)
	if (sum(tb > 1) > 0) {
		for (i in 1:sum(tb > 1)) {
			ind = which(as.character(explanatory) == names(tb)[which(tb > 1)[i]] )
			response= c(response, mean(response[ind]))
			response = response[-1 * ind]
			explanatory= explanatory[-1 * ind]
			explanatory= c(explanatory, as.numeric(names(tb)[which(tb > 1)[i]]))
		}
	}

	ind <- order(explanatory, decreasing=FALSE)
	if (sum(diff(ind) < 0) != 0) {
		explanatory <- explanatory[ind]
		response <- response[ind]
	}

	n <- length(response)
	a <- b <- eps <- rep(0, n-2)
	for (i in 2:(n-1)) {
		x <- explanatory[(i-1):(i+1)]	
		a[i-1] <- (x[3] - x[2]) / (x[3] - x[1]) 
		b[i-1] <- (x[2] - x[1]) / (x[3] - x[1]) 
		eps[i-1] <- a[i-1] * response[i-1] + b[i-1] * response[i+1] - response[i]
	}
	cSq <- 1/(a^2 + b^2 + 1)	
	list( sigmaSq=1/(n-2) * sum(cSq * eps^2), a=a, b=b, eps=eps)
}
