# if time = t[i] dann wird g im interval i-1 berechnet.
qfuncsky <- function(time,t,l,mu,psi,rho) {
	i <- interstt(time,t)	
	bivar <- bi(i,t,l,mu,psi,rho)
	aivar <- ai(i,l,mu,psi)
	top <- 4* exp(aivar*(-t[i]+time))
	bottom <- exp(aivar*(-t[i]+time))*(1+bivar)+(1-bivar)
	res <- top/(bottom)^2
	res
	}	
