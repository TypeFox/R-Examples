tau2.r.solver <- function(w, Qr){
	f <- function(tau2){
		out <- sum(sqrt(1 - w/sum(w) + tau2*(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))) - Qr*sqrt(3.14159/2)
		return(out)
	}
	f <- Vectorize(f)
	tau.upp <- Qr*sqrt(3.14159/2)/sum(sqrt(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))
	tau2.upp <- tau.upp^2
	f.low <- f(0)
	f.upp <- f(tau2.upp)
	if(f.low*f.upp > 0){
		tau2.r <- 0
	}else{
		tau2.r <- uniroot(f, interval = c(0, tau2.upp))
		tau2.r <- tau2.r$root
	}
	return(tau2.r)
}