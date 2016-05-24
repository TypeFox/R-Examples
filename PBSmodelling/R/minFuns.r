#------------------------
# Minimization Functions
#------------------------

#calcMin--------------------------------2012-12-20
calcMin <- function(pvec, func, method="nlm", trace=0, maxit=1000, reltol=1e-8, steptol=1e-6, temp=10, repN=0, ...) {
	Sfun <- function(S,pvec,func,repN) { # function of surrogate parameters
		#eval(parse(text="PBSmin$N <<- PBSmin$N + 1"))
		tget(PBSmin); PBSmin$N <- PBSmin$N + 1
		P <- restorePar(S,pvec);
		Uval <- func(P);
		if (PBSmin$N==1) {
			#eval(parse(text="PBSmin$fmin0 <<- Uval"))
			PBSmin$fmin0 <- Uval
			if (repN>0) cat("\nParameter reporting:\n")
		}
		if (repN>0 && is.element(PBSmin$N,seq(repN,maxit,repN)))
			print(paste("N=",PBSmin$N," P=c(",paste(show0(round(P,5),5),collapse=","),") F=",show0(round(Uval,5),5),sep="" ));
		tput(PBSmin)
		return(Uval); };

	Sval <- scalePar(pvec); nS <- length(Sval); junk <- gc(FALSE); 
	#eval(parse(text="PBSmin <<- list(); PBSmin$N <<- 0"))
	PBSmin <- list(); PBSmin$N <- 0; tput(PBSmin) # initialize PBSmin list
	Ftime <- proc.time()[1:3];
	if (method=="nlm") { # Non-Linear Minimization
		Fout <- nlm(f=Sfun,p=Sval,typsize=rep(1,nS), iterlim=maxit, gradtol=reltol, steptol=steptol,
		        pvec=pvec, func=func,repN=repN,...)
		Ftime <- proc.time()[1:3] - Ftime;
		Pest <- Fout$estimate; grad <- Fout$gradient; code <- Fout$code; 
		iter <- Fout$iterations; eval<- NULL; fmin <- Fout$minimum;
		mess <- switch(code,"Relative gradient is close to zero, current iterate is probably solution.",
			"Successive iterates within tolerance, current iterate is probably solution.",
			"Last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.",
			"Iteration limit exceeded.",
			"Maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.") }
	else if (method=="nlminb") { # Optimization using PORT routines
		Fout <- nlminb(start=Sval,objective=Sfun,scale=rep(1,nS),
		        control=list(trace=trace,iter.max=maxit,rel.tol=reltol,x.tol=steptol),pvec=pvec,func=func,repN=repN,...)
		Ftime <- proc.time()[1:3] - Ftime;
		Pest <- Fout$par; grad <- NULL; code <- Fout$convergence; mess <- Fout$message; 
		iter <- Fout$iterations; eval<- Fout$evaluations; fmin <- Fout$objective; }

	else { # General-purpose Optimization
		Fout <- optim(par=Sval,fn=Sfun,method=method,
		        control=list(trace=trace,maxit=maxit,reltol=reltol,temp=temp),pvec=pvec,func=func,repN=repN,...);
		Ftime <- proc.time()[1:3] - Ftime;
		omess <-c("Successful convergence.","Iteration limit maxit had been reached.",
			"Degeneracy of the Nelder-Mead simplex.","Warning from the L-BFGS-B method.",
			"Error from the L-BFGS-B method."); names(omess) <- c(0,1,10,51,52);
		Pest <- Fout$par; grad <- NULL; code <- Fout$convergence; 
		iter <- Fout$counts; eval<- NULL; fmin <- Fout$value;
		mess <- paste(omess[as.character(code)],Fout$message,sep=" ");
	}

	Pfin <- restorePar(Pest,pvec); Pmat <- cbind(Pfin,pvec[,2:4]); 
	P0 <- pvec[,1]; names(P0) <- dimnames(pvec)[[1]]; AIC <- 2*fmin + 2*nS;
	#eval(parse(text="PBSmin <<- c(PBSmin,list(start=P0, end=Pfin, surrogates=Pest, check=scalePar(Pmat), gradient=grad,
		#code=code, message=mess, iterations=iter, evaluations=eval, time=Ftime, fmin=fmin, AIC=AIC ))"))
	tget(PBSmin)
	PBSmin <- c(PBSmin,list(start=P0, end=Pfin, surrogates=Pest, check=scalePar(Pmat), gradient=grad,
		code=code, message=mess, iterations=iter, evaluations=eval, time=Ftime, fmin=fmin, AIC=AIC ))
	tput(PBSmin)
	Obag <- list(Fout=Fout, iters=iter[1], evals=PBSmin$N, cpuTime=Ftime[1], elapTime=Ftime[3], 
		fminS=PBSmin$fmin0, fminE=fmin, Pstart=P0, Pend=Pfin, AIC=AIC, message=mess);
	return(Obag) 
};
#------------------------------------------calcMin


#scalePar-------------------------------2012-12-20
scalePar <- function(pvec) { # Convert true parameters to surrogates
	Pval <- pvec[,1]; Pmin <- pvec[,2]; Pmax <- pvec[,3]; idx <- pvec[,4];
	Sval <- (Pval[idx]-Pmin[idx]) / (Pmax[idx]-Pmin[idx]);
	Sval <- pmax(Sval,0); Sval <- pmin(Sval,1);  # enforces the range
	S    <- (2/pi) * asin(sqrt(Sval)); names(S) <- dimnames(pvec)[[1]][idx];
	return(S);
}
#-----------------------------------------scalePar


#restorePar-----------------------------2012-12-20
restorePar <- function(S,pvec) { # Convert surrogates to true parameters
	Pval <- pvec[,1]; Pmin <- pvec[,2]; Pmax <- pvec[,3]; idx <- pvec[,4];
	if (sum(idx) != length(S)) stop("Warning: S & P not consistent/n");
	Pcon <- Pmin[idx] + (Pmax[idx]-Pmin[idx])*sin(pi*S/2)^2;
	P <- Pval; P[idx] <- Pcon; names(P) <- dimnames(pvec)[[1]];
	return(P); 
};
#---------------------------------------restorePar


#GT0------------------------------------2012-12-20
GT0 <- function (x, eps = 1e-04) {
    eps2 <- eps/2
    ifix <- ((x > 0) & (x < eps))
    i0 <- (x <= 0)
    y = x
    y[i0] <- eps2
    y[ifix] <- eps2 * (1 + (x[ifix]/eps)^2)
    return(y) 
}
#----------------------------------------------GT0


#===== THE END ===================================

