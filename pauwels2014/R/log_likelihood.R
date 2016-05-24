log_likelihood <-
function(simu, simu_subset, data, fit = F){
	## Likelihood
	## The fit argument is used to guide the local seach toward
	## parameters producing smooth, non exploding simulations
	## This is used to find a starting point for MH algorithm
	## in order to avoid the BFGS algorithm to reach 
	## parameter values producing badly behaved simulations
	
	if ((sum(is.na(simu)) > 0)){
		res <- -Inf
	}else{
	
		## Original model
		res <- - sum( ( (simu_subset - data)^2 / ( (0.01 + (0.2 * simu_subset)^2) ) + log( (0.01 + (0.2 * simu_subset)^2)) )[,-1] ) / 2 
		
		if(!fit){
			## Original model gives reasonable Likelihood when the simulation values are very far from original data, we penalyze this
			res <- res - sum( (simu_subset - data)^2 / (10 ) + log(10) ) / 2
		
			## We do not observe all concentrations, however, we do not want any of them to be too big
			res <- res - 1 / 2 * max( simu[,-1]^2 / 1e5 + log(1e-5) ) 

			## Penalyze large variations in concentration and abrupt discontinuities in evolution
			n <- nrow(simu)
			derivs <- simu[-c(1,2),-1] - simu[-c(n-1,n),-1]
			n1 <- nrow(derivs)
			second_derivs <- (derivs[-c(1,2),] - derivs[-c(n1-1,n1),])/2
		
			res <- res - 1/2 * (max(derivs^2) + max(second_derivs^2) )
		}
		
		if (is.na(res)){
			res <- -Inf
		}
	}
	res
}
