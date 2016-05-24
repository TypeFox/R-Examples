# Author: cns
# dlm backend

dlmodeler.filter.dlm <-
		function(yt, model, raw.result=FALSE, logLik=FALSE, filter=TRUE)
{
	if(!require('dlm')) stop("required package could not be found: dlm")
	Ttvar <- length(dim(model$Tt))==3
	Rtvar <- length(dim(model$Rt))==3
	Qtvar <- length(dim(model$Qt))==3
	Ztvar <- length(dim(model$Zt))==3
	Htvar <- length(dim(model$Ht))==3
	if( Ttvar|Rtvar|Qtvar|Htvar ) stop("This kind of time varying model is currently unsupported")
	# TODO handle more time varying cases
	Ht <- model$Rt %*% model$Qt %*% t(model$Rt)
	if( Ztvar ) {
		n1 <- NROW(model$Zt)
		n2 <- NCOL(model$Zt)
		n3 <- dim(model$Zt)[3]
		FF <- matrix(0,n1,n2)
		JFF <- matrix(1:(n1*n2),n1,n2)
		X <- matrix(NA,n3,n1*n2)
		for( i in 1:n1 ) for( j in 1:n2 ) X[,i+(j-1)*n1] <- model$Zt[i,j,]
		mdlm <- dlm::dlm(
				m0=model$a0,
				C0=model$P0+1e7*model$P0inf,
				FF=FF,JFF=JFF,X=X,
				V=model$Ht,
				GG=model$Tt,
				W=Ht)
	} else {
		mdlm <- dlm::dlm(
				m0=model$a0,
				C0=model$P0+1e7*model$P0inf,
				FF=model$Zt,
				V=model$Ht,
				GG=model$Tt,
				W=Ht)
	}
	# only filter if needed
	if( filter ) {
		res <- dlm::dlmFilter(t(yt),mdlm,simplify=!raw.result)
		res.Pt <- dlm::dlmSvd2var(res$U.R,res$D.R)
		Pt <- array(NA,dim=c(NROW(Ht),NCOL(Ht),length(res.Pt)))
		for( i in 1:length(res.Pt) ) Pt[,,i] <- res.Pt[[i]]
    f <- t(res$f)
    at <- t(res$a)
	} else {
	  f <- NA
	  at <- NA
	  Pt <- NA
	  res <- NA
	}
  # only compute log likelihood if needed
	# note: dlm currently does not compute the log likelihood when filtering
  # so it needs to be computed separately
	if( logLik ) logLik <- dlm::dlmLL(yt,mdlm) else logLik <- NA
	if( raw.result ) raw.res <- res else raw.res <- NA
	
	return(list(backend='dlm',
					f=f,
					at=at,
					Pt=Pt,
					logLik=logLik,
					d=0,
					raw.result=raw.res))
}



dlmodeler.smooth.dlm <-
		function(filt, raw.result=FALSE)
{
	if(!require('dlm')) stop("required package could not be found: dlm")
	
	res <- dlm::dlmSmooth.dlmFiltered(filt$raw.result)
	
	if( raw.result ) raw.res <- res else raw.res <- NA
	res.Pt <- dlm::dlmSvd2var(res$U.S,res$D.S)
	Pt <- array(NA,dim=c(NCOL(res$D.S),NCOL(res$D.S),length(res.Pt)))
	for( i in 1:length(res.Pt) ) Pt[,,i] <- res.Pt[[i]]
	
	return(list(backend='dlm',
					at=t(res$s),
					Pt=Pt,
					raw.result=raw.res))
}


