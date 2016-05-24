hazards  <-  function(dist="Weibull", t, parms){

if(dist=="Weibull")	 	haz <- (parms[1]^parms[2])*parms[2]*t^(parms[2]-1) 
else if(dist=="loglogistic")	haz <- parms[1]*parms[2]*t^(parms[2]-1)/(1+parms[1]*t^parms[2])
else if(dist=="Gompertz") haz <- parms[1]*exp(parms[2]*t)
else if(dist=="lognormal") haz <- dnorm((log(t)-parms[1])/parms[2])/(parms[2]*t*(1-pnorm((log(t)-parms[1])/parms[2])))
else if(dist=="gamma") haz <- dgamma(t,shape=parms[2], scale=1/parms[1])/(1-pgamma(t,shape=parms[2], scale=1/parms[1]))
else stop("Unrecognized baseline distribution")
 haz
 }