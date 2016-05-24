# Author: cns
# FKF backend

dlmodeler.filter.FKF <-
		function(yt, model, raw.result=FALSE, logLik=TRUE, filter=TRUE)
{
	if(!require('FKF')) stop("required package could not be found: FKF")
	rqr.fun <- function(r,q) r %*% q %*% t(r)
	RQt <- dlmodeler.timevar.fun(model$Rt, model$Qt, rqr.fun)
	m <- NROW(model$Tt)
	d <- NROW(model$Zt)
	res <- FKF::fkf(
			a0=as.vector(model$a0),
			P0=model$P0+1e7*model$P0inf,
			dt=matrix(0,m,1),ct=matrix(0,d,1),
			Tt=model$Tt,
			Zt=model$Zt,
			HHt=RQt,
			GGt=model$Ht,
			yt)
	if( raw.result ) raw.res <- res else raw.res <- NA
	if( length(dim(model$Zt))==2 ) {
		# 1-step ahead prediction when observation matrix is not time-varying
		f <- model$Zt %*% res$at
	} else {
		# 1-step ahead prediction when observation matrix is time-varying
		f <- matrix(NA,NROW(model$Zt),NCOL(yt))
		for( i in 1:NCOL(yt) ) f[,i] <- model$Zt[,,i] %*% res$at[,i]
	}
	return(list(backend='FKF',
					f=f,
					at=res$at,
					Pt=(res$Pt),
					logLik=res$logLik,
					d=0,
					raw.result=raw.res))
}



dlmodeler.smooth.FKF <-
		function(filt, raw.result=FALSE)
{
	stop("smoothing is not supported for backend: FKF");
}

