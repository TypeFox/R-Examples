min.lnorm <-
function(p, x, delta){
	meanlog <- p[1]
	sdlog <- p[2]
	LL <- sum(delta*log(dlnorm(x, meanlog=meanlog,sdlog=sdlog))) + sum((1-delta)*log(1-plnorm(x, meanlog=meanlog,sdlog=sdlog)))
	return(-LL)
	}

