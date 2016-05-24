Update_aE <-
function(last.params){
		aE_prime <- last.params$aE+rnorm(length(last.params$aE),0,last.params$aE_stp)
		prior_prob_alphaE_prime <- Prior_prob_alphaE(aE_prime) 
		new.params <- last.params
		if(prior_prob_alphaE_prime != -Inf) {
			covariance_prime <- Covariance(last.params$a0,last.params$aD,aE_prime,last.params$a2,last.params$D,last.params$E,last.params$delta)
			if(is.positive.definite(covariance_prime)){
				LnL_thetas_vec_prime <- Likelihood_thetas(last.params$thetas,covariance_prime)
					if(exp((prior_prob_alphaE_prime+sum(LnL_thetas_vec_prime)) - (last.params$prior_prob_alphaE+sum(last.params$LnL_thetas_vec))) >= runif(1)){
						new.params$aE <- aE_prime
						new.params$covariance <- covariance_prime
						new.params$prior_prob_alphaE <- prior_prob_alphaE_prime						
						new.params$LnL_thetas_vec <- LnL_thetas_vec_prime
						new.params$aE_accept <- new.params$aE_accept + 1 				
					}
			}	
		}
        new.params$aE_moves <- new.params$aE_moves + 1		
        return(new.params)
	}
