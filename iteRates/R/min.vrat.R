min.vrat <-
function(p, x, delta){
	alpha <- p[1]
	beta <- p[2]
	LL <- sum(delta*(log(alpha*beta)-(1+alpha)*log(1+beta*x))) + sum((1-delta)*log((1+beta*x)^(-alpha)))
	return(-LL)
	}

