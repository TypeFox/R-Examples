censor.exp.x <-
function(delta,x,min.branch){
	Fn <- sum(delta)
	x[x==0] <- min.branch
	Tn <- sum(x)
	est <- Fn/Tn
	LL <- Fn*(log(Fn)-log(Tn)-1)
	a <- data.frame(t(c(est,NA)),LL)
	return(a)
	}

