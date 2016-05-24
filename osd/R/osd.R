osd <- function(D, k=NULL, y.profile=NULL, res.method = c("ica.osd","osd","icr","mcr"), comp.coef = 2, noise.thr=0.01)
{

	res.method <- match.arg(res.method, c("ica.osd","osd","icr","mcr"))

	D.int <- D
	D[is.na(D)] <- 0
	win.maxs <- apply(D,2,function(x){max(x)})
	D[,win.maxs<max(D)*noise.thr] <- 0
				
	D.comp <- D^(1/comp.coef)
	D.comp[is.na(D.comp)] <- 0
	
	if(res.method=="icr")
	{
		if(is.null(k)) stop("k is missing. Please, specify the number of factors to deconvolve.")
		ICA.res <- ica.resolution(D.comp, k)
		Cn <- chrom.isoreg(ICA.res$S)^comp.coef
		
		osd.res <- list(C=Cn,S=matrix(0,ncol=ncol(Cn),nrow=ncol(D)))
		
		nonreg.res <- nonnegreg(Cn, osd.res$S, D)
		osd.res$S <- nonreg.res$S
		osd.res$S[normalize(osd.res$S)<0.005] <- 0
		osd.res$C <- nonreg.res$C
	}
	if(res.method=="mcr")
	{
		if(is.null(k)) stop("k is missing. Please, specify the number of factors to deconvolve.")
		MCR.res <- mcr.resolution(D.comp, k)

		Cn <- chrom.isoreg(MCR.res$C)^comp.coef
		osd.res <- list(C=Cn,S=MCR.res$S,nrow=ncol(D))
		
		if(comp.coef!=1)
		{
			nonreg.res <- nonnegreg(osd.res$C, osd.res$S, D)
			osd.res$S <- nonreg.res$S
			osd.res$C <- nonreg.res$C	
		} 
		osd.res$S[normalize(osd.res$S)<0.005] <- 0
	}
	if(res.method=="ica.osd")
	{
		if(is.null(k)) stop("k is missing. Please, specify the number of factors to deconvolve.")
		ICA.res <- ica.resolution(D.comp, k)
		Cn <- chrom.isoreg(ICA.res$S)^comp.coef
				
		osd.res <- getS_wOSD(Cn, matrix(0,ncol=ncol(Cn),nrow=ncol(D)), D)	
		osd.res$S[normalize(osd.res$S)<0.005] <- 0	
	}
	if(res.method=="osd")
	{
		if(is.null(y.profile)) stop("y.profile is missing. Please, specify the compound elution profile model.")
		Sn <- getS.OSD(D=D, mod.c=y.profile, beta=comp.coef)	
		Sn[normalize(Sn)<0.005] <- 0	
		
		osd.res <- list(C=NULL,S=Sn)
		
	}	
	
	osd.ret <- list(data = D.int, C = osd.res$C, S =osd.res$S, k=k, y.profile=y.profile, res.method = res.method, comp.coef = comp.coef, noise.thr=noise.thr)
	
	#osd.ret <- new(Class="OSDres", data = D.int, C = osd.res$C, S =osd.res$S, k=as.integer(k), y.profile=as.vector(y.profile), res.method = res.method, comp.coef = comp.coef, noise.thr=noise.thr)
	osd.ret
}
