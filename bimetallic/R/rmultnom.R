rmultnom <-
function(N, probs){
	ret = probs * 0
	tab = table(findInterval(runif(N), cumsum(probs)))
	ret[as.numeric(names(tab)) + 1] = tab
	ret
	}

