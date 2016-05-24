Prior_prob_alphaE <-
function(aE){
		sum( dexp(aE,log=TRUE) )
	}
