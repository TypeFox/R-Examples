min.weib <-
function(p, Fn, n, x, delta){
	theta <- p[1]
	k <- p[2]
	LL <- Fn*log(k) + k*Fn*log(theta) - sum((x*theta)^k) + (k-1)*sum(delta*log(x))
	return(-LL)
	}

