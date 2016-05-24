Update_aD <-
function(last.params){
		aD_prime <- last.params$aD+rnorm(1,0,last.params$aD_stp)
		prior_prob_alphaD_prime <- Prior_prob_alphaD(aD_prime) 
		new.params <- last.params
		if(prior_prob_alphaD_prime != -Inf) {											
			covariance_prime <- Covariance(last.params$a0,aD_prime,last.params$aE,last.params$a2,last.params$D,last.params$E,last.params$delta)
			if(is.positive.definite(covariance_prime)){
				LnL_thetas_vec_prime <- Likelihood_thetas(last.params$thetas,covariance_prime) 		
					if(exp((prior_prob_alphaD_prime+sum(LnL_thetas_vec_prime)) - (last.params$prior_prob_alphaD+sum(last.params$LnL_thetas_vec))) >= runif(1)){		
						new.params$aD <- aD_prime
						new.params$covariance <- covariance_prime
						new.params$prior_prob_alphaD <- prior_prob_alphaD_prime						
						new.params$LnL_thetas_vec <- LnL_thetas_vec_prime
						new.params$aD_accept <- new.params$aD_accept + 1 				
					}
			}		
		}
        new.params$aD_moves <- new.params$aD_moves + 1		
        return(new.params)
	}
