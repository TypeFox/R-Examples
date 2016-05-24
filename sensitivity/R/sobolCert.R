# Certified Sobol SA
# Author: Alexandre Janon <alexandre.janon at imag.fr>
# see also src/as.cpp and src/rglueas.cpp
# Reference: Janon, A., Nodet M., Prieur C. (2011) Uncertainties assessment in global sensitivity indices estimation from metamodels.

genDesignPoints=function(p,N,X1,X2) {
	.C("Rglue_set_as", as.integer(p), as.integer(N), as.double(X1), as.double(X2))
	M=rep(0,p*N*(p+1))
	M=.C("Rglue_get_design",A=as.double(M))$A
	return(matrix(M,ncol=p,nrow=N*(p+1),byrow=T))
}


sobolCert=function(model=NULL, X1=NULL, X2=NULL, nboot=300, conf=0.95, lambda0=0, h=0) {
	x=list(call=match.call())
	class(x)="sobolCert"
	.C("Rglue_set_lambda0_r",as.double(lambda0),as.double(h))
	if(!is.null(X1)) {
		p=ncol(X1)
		N=nrow(X1)
		if(ncol(X2)!=p || nrow(X2)!=N) {
			stop("The samples X1 and X2 must have the same dimensions")
		}
		x$X=genDesignPoints(p,N,t(X1),t(X2))
		for(i in 1:nrow(x$X)) {
			param=x$X[i,]
			m=model(param)
			x$y[i]=m$out
			x$err[i]=m$err
		}
		.C("Rglue_fill_vf",as.double(x$y),as.double(x$err))
	} else {
		p=0
		N=0
		out=.C("Rglue_getPN",p=as.integer(p),N=as.integer(N))
		p=out$p
		N=out$N
	}
	minn=rep(0,p)
	maxn=rep(0,p)
	valmin=rep(0,p)
	valmax=rep(0,p)

	if(lambda0<1e-5) {
		sob=.C("Rglue_sobol_regr_bc", as.integer(nboot), Smin=as.double(minn), Smax=as.double(maxn), valmin=as.double(valmin), valmax=as.double(valmax),risque=as.double(1-conf))
	} else {
		x$lambda0=lambda0
		x$h=h
		sob=.C("Rglue_sobol_BFGS", as.integer(nboot), Smin=as.double(minn), Smax=as.double(maxn), valmin=as.double(valmin), valmax=as.double(valmax),risque=as.double(1-conf))
		x$penalty=sobolCert.penalty()
	}
	x$S=sob
	return(x)
}

sobolCert.penalty=function(lambda=0, h=0) {
	p=0
	N=0
	out=.C("Rglue_getPN",p=as.integer(p),N=as.integer(N))
	p=out$p
	N=out$N
	valopti=rep(0,p)
	pmin=rep(0,p)
	pmax=rep(0,p)
	z=rep(0,2*N)

	if(lambda>0) {
		.C("Rglue_set_lambda0_r",as.double(lambda),as.double(h))
	}

	for(i in 0:(p-1)) {
		print(i)
		penal=0
		sob=.C("Rglue_sobol_BFGS_replica", valopti=as.double(valopti), as.integer(i), as.integer(1), as.integer(-1), z=as.double(z),p=as.double(penal))
		pmin[i+1]=sob$p
	#	sob=.C("Rglue_sobol_BFGS_replica", valopti=as.double(valopti), as.integer(i), as.integer(-1), as.integer(-1), z=as.double(z),p=as.double(penal))
	#	pmax[i+1]=sob$p
	}
#	out=list(pmin=pmin, pmax=pmax)
	return(pmin)
}


print.sobolCert=function(x, ...) {
   cat("\nCall:\n", deparse(x$call), "\n", sep = "")
	S=x$S
	Sdf=data.frame(ptmin=S$valmin, ptmax=S$valmax, smin=S$Smin, smax=S$Smax)
	colnames(Sdf)=c("low. bnd.", "upp. bnd.", "min. c.i.", "max. c.i.")
	print(Sdf)
}

