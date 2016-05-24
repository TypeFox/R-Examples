p <- function(i,time,t,l,mu,psi,rho) {
	if (t[i]>time) {res = 0} else{
		bivar <- bi(i,t,l,mu,psi,rho)
		aivar <- ai(i,l,mu,psi)
		temp<-exp(aivar*(-t[i]+time))*(1+bivar)
		res <- (l[i]+mu[i]+psi[i]-aivar*(temp-(1-bivar))/(temp+(1-bivar)))
		res<- res/(2*l[i])
	}
	res
	}
