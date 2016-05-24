Covariance <-
function(a0,aD,aE,a2,GeoDist,EcoDist,delta) {
        weighted.dist.squared <- aE[1] * (EcoDist[[1]])^2
        if (length(EcoDist)>1) {
            for (k in 2:length(EcoDist)) {
                weighted.dist.squared <- weighted.dist.squared + aE[k] * (EcoDist[[k]])^2
            }
        }
		covariance <- (1/a0)*exp((-sqrt(aD*GeoDist^2 + weighted.dist.squared)^a2)-delta)
		diag(covariance) <- (1/a0)
		return(covariance)
	}