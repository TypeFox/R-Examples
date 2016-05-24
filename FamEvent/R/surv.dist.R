surv.dist <- function(base.dist, t, parms, xbeta, alpha, res){
	if(base.dist=="Weibull") exp(-(parms[1]*t)^parms[2]*exp(xbeta+alpha))-res
	else if(base.dist=="loglogistic") 1/(1+parms[1]*t^parms[2])*exp(xbeta+alpha) -res
	else if(base.dist=="Gompertz") exp(-(parms[1]*(exp(parms[2]*t)-1)/parms[2])*exp(xbeta+alpha)) -res
	else if(base.dist=="lognormal"){
		x <- (log(t)-parms[1])/parms[2]
		H <- -log(1-pnorm(x))
		exp(-H*exp(xbeta+alpha))-res
	}
	else if(base.dist=="gamma"){
		H <- -log(1-pgamma(t,shape=parms[1], scale=parms[2] ))
		exp(-H*exp(xbeta+alpha)) - res
	}
	else stop("Unrecognized baseline distribution")
} 
  	