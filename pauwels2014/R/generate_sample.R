generate_sample <-
function(theta, knobj, N = 500, step = 1, verbose = F){
	## Simple Metropolis hasting implementation
	## theta is an initial value
	## Proposal is gaussian with variance 1
	## with possible single coordinate updates
	## Variance is decreased when the rejection rate is too high
	
	#print(seed)
	#set.seed(seed)
	
	thetas <- matrix(0,N, length(theta))
	thetas[1,] <- theta
	lls <- rep(0,N)
	lls[1] <- eval_log_like_knobj(theta, knobj, fit = T)
	n_accept <- 1
	
	for (i in 2:(N)){
		
		## Gaussian proposal
		theta2 <- theta + step * rnorm(length(theta))
		ll <- eval_log_like_knobj(theta2, knobj, fit = T)
		u <- runif(min=0,max=1,n=1)
		
		if (u < exp(ll - lls[i-1])){
			thetas[i,] <- theta2
			theta <- theta2
			lls[i] <- ll
			n_accept <- n_accept + 1
		}else{
			thetas[i,] <- theta
			lls[i] <- lls[i-1]
		}
		
		## Single coordinate proposal
		temp <- sample(1:length(theta),size = 1)
		theta2 <- theta
		theta2[temp] <- theta[temp] + rnorm(1)
		ll <- eval_log_like_knobj(theta2, knobj, fit = T)
		u <- runif(min=0,max=1,n=1)
		
		if (u < exp(ll - lls[i-1])){
			thetas[i,] <- theta2
			theta <- theta2
			lls[i] <- ll
		}else{
			thetas[i,] <- theta
			lls[i] <- lls[i-1]
		}
		
		
		if((i%%100 == 0)&(i < (N/2))){
			## Decrease proposal variance if rejection rate is too high
			if((n_accept/100) < 0.0099){
				n_accept = 1
				step <- step/2
				print(paste("Divide step:", step))
			}else{
				n_accept = 1
			}
		}
		if((i%%100 == 0)&verbose){
			print(i)
		}
	}
	
	n_accept / N
	dimnames(thetas)[[2]] <- knobj$global_parameters$param_names
	thetas
}
