Prior_prob_mu <-
function(parameter_mu,beta){
		dnorm(parameter_mu[1,],0,sd=sqrt(1/beta),log=TRUE)
	}
