g <-
function(time,t,l,mu,rho) {
	i <- inter(time,t)	
	ci <- const(i,t,l,mu,rho)	
	zah <- exp(-(l[i]-mu[i]) * (t[i]+time))
	a<-(ci-1) * exp(-(l[i]-mu[i]) * t[i])
	b<- (mu[i] - ci * l[i]) * exp(-(l[i]-mu[i]) * time) 
	nen <- l[i] * a +  b
	res <- zah/(nen)^2
	res
	}

