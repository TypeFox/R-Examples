Update_a2 <-
function(last.params){
		a2_prime <- last.params$a2+rnorm(1,0,last.params$a2_stp)
		prior_prob_alpha2_prime <- Prior_prob_alpha2(a2_prime) 
		new.params <- last.params
		if(prior_prob_alpha2_prime != -Inf) {
			covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$aE,a2_prime,last.params$D,last.params$E,last.params$delta)
			if(is.positive.definite(covariance_prime)){
				LnL_thetas_vec_prime <- Likelihood_thetas(last.params$thetas,covariance_prime)
					if(exp((prior_prob_alpha2_prime+sum(LnL_thetas_vec_prime)) - (last.params$prior_prob_alpha2+sum(last.params$LnL_thetas_vec))) >= runif(1)){
						new.params$a2 <- a2_prime
						new.params$covariance <- covariance_prime	
						new.params$prior_prob_alpha2 <- prior_prob_alpha2_prime
						new.params$LnL_thetas_vec <- LnL_thetas_vec_prime
						new.params$a2_accept <- new.params$a2_accept + 1 				
					}
			}		
		}
	new.params$a2_moves <- new.params$a2_moves + 1	
	return(new.params)
	}
