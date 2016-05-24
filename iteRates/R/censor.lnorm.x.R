censor.lnorm.x <-
function(delta,x,min.branch){
	x[x==0] <- min.branch
	est <- optim(c(0,1), min.lnorm, delta=delta, x=x,lower=c(-Inf,.Machine$double.xmin),method="L-BFGS-B")
	par <- est$par
	LL <- -est$value
	
	a <- data.frame(t(par),LL)
	return(a)
	}

