Update_beta <-
function(last.params){
		new.params <- last.params
		new.params$beta <- rgamma(1,shape=(last.params$loci/2+0.001),rate=(0.001+sum(((last.params$mu[1,])^2)/2)))  
		new.params$prior_prob_beta <- Prior_prob_beta(new.params$beta)
		new.params$prior_prob_mu <- Prior_prob_mu(last.params$mu,new.params$beta)
		new.params$beta_moves <- new.params$beta_moves+1
		return(new.params)
	}

