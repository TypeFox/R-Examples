fhatgg <-
function(x,theta,mu,sigma){
	dn=theta*dnorm(x,mu,sigma)
	f=sum(dn)
	fg=sum(dn*(mu-x)/sigma^2)
	fgg=sum(dn*((mu-x)^2/sigma^2-1)/sigma^2)
	return(c(f,fg,fgg))
}
