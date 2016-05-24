Update_a0 <-
function(last.params){
		new.params <- last.params
			gamma_rate_part <- a0_gibbs_rate(last.params$thetas,last.params$covariance,last.params$a0)
		new.params$a0 <- rgamma(1,shape=(1+(last.params$loci*last.params$k)/2),rate=(1+gamma_rate_part))	
		new.params$covariance <- Covariance(new.params$a0,last.params$aD,last.params$aE,last.params$a2,last.params$D,last.params$E,last.params$delta)
		new.params$prior_prob_alpha0 <- Prior_prob_alpha0(new.params$a0)
		new.params$LnL_thetas_vec <- Likelihood_thetas(last.params$thetas,new.params$covariance)
		new.params$a0_moves <- last.params$a0_moves + 1
		return(new.params)
	}
