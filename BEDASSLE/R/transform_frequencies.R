transform_frequencies <-
function(thetas,mu){
		allele.frequencies <- 1/(1+exp(-(thetas+mu)))
		return(allele.frequencies)
	}
