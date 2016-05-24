tau2.m.solver <- function(w, Qm){
	f <- function(tau2){
		out <- sum(sqrt(1 + w*tau2)) - Qm*sqrt(3.14159/2)
		return(out)
	}
	f <- Vectorize(f)
	n <- length(w)
	tau2.upp <- sum(1/w)*(Qm^2/n*2/3.14159 - 1)
	tau2.upp <- max(c(tau2.upp, 0.01))
	f.low <- f(0)
	f.upp <- f(tau2.upp)
	if(f.low*f.upp > 0){
		tau2.m <- 0
	}else{
		tau2.m <- uniroot(f, interval = c(0, tau2.upp))
		tau2.m <- tau2.m$root
	}
	return(tau2.m)
}