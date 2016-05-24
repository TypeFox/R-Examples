censor.weib.x <-
function(delta,x,min.branch){
	x[x==0] <- min.branch
	Fn <- sum(delta)
	n <- length(x) 

	est <- optim(c(1,1), min.weib, Fn=Fn, delta=delta, x=x, n=n)
	par <- est$par
	LL <- -est$value
	
	a <- data.frame(t(par),LL)
	return(a)
	}

