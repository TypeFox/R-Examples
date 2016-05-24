# Author: cns
# KFAS backend

dlmodeler.filter.KFAS <-
		function(yt, model, raw.result=FALSE, logLik=TRUE, filter=TRUE)
{
	# KFAS is a dependency
	# if(!require('KFAS')) stop("required package could not be found: KFAS")
	kfas.model <- KFAS::SSModel( t(yt) ~ -1 + SSMcustom(
	  Z=model$Zt, # observation
	  T=model$Tt, # transition
	  R=model$Rt, # state disturbance selection
	  Q=model$Qt, # state disturbance covariance
	  a1=model$a0, # initial state
    P1=model$P0, # initial state covariance
    P1inf=model$P0inf, # diffuse part of P1
    n=NROW(yt)),
    H=model$Ht # observation disturbance
  )
  
	# only filter if needed
  if( filter ) {
  	res <- KFAS::KFS(kfas.model,smoothing="none") # just filtering
  	if( length(dim(model$Zt))==2 ) {
  	  # 1-step ahead prediction when observation matrix is not time-varying
  	  f <- kfas.model$Z[,,1] %*% t(res$a)
  	} else {
  	  # 1-step ahead prediction when observation matrix is time-varying
  	  f <- matrix(NA, NROW(res$model$Z), NCOL(yt))
  	  for( i in 1:NCOL(yt) ) f[,i] <- res$model$Z[,,i] %*% res$a[i,]
  	}
    at <- t(res$a)
    Pt <- res$P
    d <- res$d
  } else {
    f <- NA
    at <- NA
    Pt <- NA
    res <- NA
    d <- NA
  }
	# only compute log likelihood if needed
	# note: apparently KFAS does not compute the log likelihood when filtering
	# so it needs to be computed separately
  if( logLik ) logLik <- logLik(kfas.model) else logLik <- NA
	if( raw.result ) raw.res <- res else raw.res <- NA
	
	return(list(backend='KFAS',
					f=f,
					at=at,
					Pt=Pt,
					logLik=logLik,
					d=d,
					raw.result=raw.res))
}



dlmodeler.smooth.KFAS <-
		function(filt, raw.result=FALSE)
{
	# KFAS is a dependency
	# if(!require('KFAS')) stop("required package could not be found: KFAS")
	
	res <- KFAS::KFS(filt$raw.result$model,smoothing="state") # same function, now smoothing
	
	if( raw.result ) raw.res <- res else raw.res <- NA
	
	return(list(backend='KFAS',
					at=t(res$alphahat),
					Pt=res$V,
					raw.result=raw.res))
}


