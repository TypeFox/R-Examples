hppSimulate <-
function(lambda, maxVal) {
	nEvents = rpois(1, lambda*maxVal)
	hppEvents = runif(nEvents, min=0, max=maxVal)
	return(hppEvents)
}

