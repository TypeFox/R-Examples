loglik.anoint <- function(fit){
	if(any(class(fit)=="coxph")){
		fit$loglik[2]
	}	
	else{
		as.numeric(logLik(fit))
	}
}

LRT.onebyone <- function(OBO,OBO.null){

	OBO.LRTs <- sapply(OBO,loglik.anoint)
	OBO.NULL.LRTs <- sapply(OBO.null,loglik.anoint)

-2*(OBO.NULL.LRTs-OBO.LRTs)
}

LRT.uim <- function(uim,uim.null){

	uim.lrt <- loglik.anoint(uim)
	uim.null.lrt <- loglik.anoint(uim.null)

-2*(uim.null.lrt-uim.lrt)
}
