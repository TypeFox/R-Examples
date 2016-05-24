optimise.copula.mi <-
function(mis,distribution,interval=c(-10,5),npoints){
	
	#quadratic loss function
	loss.fun = function(c,mi,...){
		(mi-generate.patchwork.copula(c=exp(c),...,returnmi=TRUE)$mi)^2
	}
	
	cvals = sapply(mis,function(mi){optimize(loss.fun,interval=interval,p=distribution,mi=mi,bins=ncol(distribution),npoints=npoints)$minimum})
	return(exp(cvals))
}
