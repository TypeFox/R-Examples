q2 <-
function(i,time,t,l,mu,rho) {
	if (t[i]>time) {res = 0} else{
	ci <- const(i,t,l,mu,rho)
	a<-(ci-1) * exp(-(l[i]-mu[i]) * t[i])
	b<- (mu[i] - ci * l[i]) * exp(-(l[i]-mu[i]) * time) 
	res <- (mu[i] * a +b ) / (l[i] * a + b )
	}
	res
	}

