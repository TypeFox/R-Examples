censor.vrat.x <-
function(delta,x,min.branch){
	x[x==0] <- min.branch
	est <- optim(c(1,1),min.vrat, NULL, delta=delta, x=x,method="L-BFGS-B",lower=c(1,1e-10))
	par <- est$par
	LL <- -est$value
	
	a <- data.frame(t(par),LL)
	return(a)
	}

