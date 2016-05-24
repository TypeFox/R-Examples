# if time = t[i] dann wird g im interval i-1 berechnet.
qfuncskylog <- function(time,t,l,mu,psi,rho) {
	i <- interstt(time,t)	
	bivar <- bi(i,t,l,mu,psi,rho)
	aivar <- ai(i,l,mu,psi)
	top <- log(4) + (aivar*(-t[i]+time))
	#bottom <- log(exp(aivar*(-t[i]+time))*(1+bivar)+(1-bivar))
	bottom <- aivar*(-t[i]+time)+log((1+bivar)*(1+(1-bivar)/(exp(aivar*(-t[i]+time))*(1+bivar))))
	res <- top-2*(bottom)
	res
	}		
