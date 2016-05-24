Prior_prob_beta <-
function(beta){
		dgamma(beta,shape=0.001,rate=0.001,log=TRUE)
	}
