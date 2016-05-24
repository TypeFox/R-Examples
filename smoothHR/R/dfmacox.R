dfmacox <- function(time, time2=NULL, status, nl.predictors, other.predictors, smoother, method, mindf=NULL, maxdf=NULL, ntimes=NULL, data) {
	options(warn=-1)
	ctype <- "FALSE"
	if ( missing(data) ) stop("The argument data is missing")
	if ( missing(time) ) stop("The argument time is missing")
	if ( missing(status) ) stop("The argument status is missing")
	if ( missing(nl.predictors) ) stop("The argument 'nl.predictors' is missing")
	if ( !missing(time2) ) ctype <- "TRUE"
	#if ( missing(mcaic) ) mcaic <- "FALSE"
	if ( missing(smoother) ) smoother <- "ns"
	if (smoother != "ns" & smoother != "pspline") stop("argument 'smoother' must be 'ns' or 'pspline'")
	if ( missing(method) ) method <- "AIC"
	if ( missing(ntimes) ) ntimes <- 5
	if ( method != "AIC" & method != "AICc" & method != "BIC") stop("The argument 'method' is not valid")
	if ( ntimes < 2 | ntimes > 10 ) stop("The argument 'ntimes' must be between 2 and 10")
	p0 <- match(names(data), time, nomatch=0)
	p1 <- which(p0 == 1)
	ntime <- data[,p1]
	if (sum(p0) == 0) stop("variable defined in argument 'time' is not in the dataset 'data'")
	p2 <- match(names(data), status, nomatch=0)
	p3 <- which(p2 == 1)
	nstatus <- data[,p3]
	if (sum(p2) == 0) stop("variable defined in argument 'status' is not in the dataset 'data'")
	if (ctype == TRUE) {
		p2 <- match(names(data), time2, nomatch=0)
		if (sum(p2) == 0) stop("variable defined in argument 'time2' is not in the dataset 'data'")
		p3 <- which(p2 == 1)
		ntime2 <- data[,p3]
	}
	if ( !missing(other.predictors) ) opred <- paste(other.predictors, collapse="+")
	npred0 <- paste("pspline(", nl.predictors, ", df=0, caic=TRUE)", collapse="+")
	if ( !missing(other.predictors) ) pred <- paste(npred0, opred, sep="+")
	else pred <- npred0
	nnl <- length(nl.predictors)
	if (!missing(maxdf) & length(maxdf)!=nnl) stop("The argument 'maxdf' must be a vector with the same length as the number of nonlinear variables")
	if (missing(maxdf) & smoother == "ns") maxdf <- rep(15,nnl)
	if (missing(mindf) & smoother == "ns") mindf <- rep(1,nnl)
	if (missing(mindf) & smoother == "pspline") mindf <- rep(1.5,nnl)
	for (i in 1:nnl) {
		p2 <- match(names(data), nl.predictors[i], nomatch=0)
		if (sum(p2)==0) stop("Check variables in argument 'nl.predictors'")
	}
	if ( !missing(other.predictors) ) {
		nop <- length(other.predictors)
		for (i in 1:nop) {
			p2 <- match(names(data), other.predictors[i], nomatch=0)
			if (sum(p2) == 0) stop("Check variables in argument 'other.predictors'")
		}
	}
	if (dim(table(nl.predictors)) < nnl) stop("Check variables in argument nl.predictors")
	if ( !missing(other.predictors) ) {
		nop <- length(other.predictors)
		if (dim( table(other.predictors) ) < nop) stop("Check variables in argument other.predictors")
		all.predictors <- c(other.predictors, nl.predictors)
		if (dim( table(all.predictors) ) < nop+nnl) stop("Check variables in argument's nl.predictors and other.predictors")
	}
	if (smoother == "pspline") {
		if (nnl > 4) stop("The maximum number of nonlinear predictors is 4")
		if (nnl == 3) cat("This may take a few seconds...\n")
		if (nnl > 3) cat("This may take a few minutes...\n")
		nnl1 <- seq(1:nnl)
		ndf <- paste("df[", nnl1, "]", sep="")
		if (ctype == "TRUE") {
			p4 <- match(names(data), time2, nomatch=0)
			p5 <- which(p4 == 1)
			ntime2 <- data[,p5]
		}
		if (ctype == "TRUE") {
			covar <- as.formula( paste(" Surv(ntime,ntime2,nstatus)~ ", pred) )
			fit <- coxph(covar, data=data, x=TRUE)
		} else {
			covar <- as.formula( paste(" Surv(ntime,nstatus)~ ", pred) )
			fit <- coxph(covar, data=data, x=TRUE)
		}
		if (nnl == 1) {
			if ( is.null(time2) ) ndf1 <- dfpsuniv(time=time, status=status, nl.predictors=nl.predictors, other.predictors=other.predictors, method=method, data=data)
			else ndf1 <- dfpsuniv(time=time, time2=time2, status=status, nl.predictors=nl.predictors, other.predictors=other.predictors, method=method, data=data)
		} else {
			#Step1
			df0 <- c( rep(1.5, nnl) )
			if ( !missing(mindf) ) df0 <- mindf
			df1 <- fit$df[1:nnl]
			if ( !missing(maxdf) ) df1 <- maxdf
			df2 <- cbind(df0, (df0+df1)/2, df1)
			for (i in 1:ntimes) {
				aaa <- lapply(apply(df2, 1, function(z) {list(c(z[1],z[2],z[3]) );}), function(y) {unlist(y);})
				out <- do.call(expand.grid, aaa)
				myaic <- rep(100000, dim(out)[1])
				for (k in 1:dim(out)[1]) {
					df <- out[k,]
					npred1 <- paste("pspline(", nl.predictors, ",df=", df, ")", collapse="+")
					if ( !missing(other.predictors) ) pred1 <- paste(npred1, opred, sep="+")
					else pred1 <- npred1
					if (ctype == "TRUE") {
						covar <- as.formula( paste(" Surv(ntime,ntime2,nstatus)~ ", pred1) )
						try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
					} else {
						covar <- as.formula( paste(" Surv(ntime,nstatus)~ ", pred1) )
						try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
					}
					if (method == "AIC") try(myaic[k] <- -2*fit1$loglik[2]+2*sum(fit1$df), TRUE)
					else if (method == "AICc") try(myaic[k] <- -2*fit1$loglik[2]+2*( fit1$nevent*(sum(fit1$df)+1)/(fit1$nevent-sum(fit1$df)-2) ), TRUE)
					else if (method == "BIC") try(myaic[k] <- -2*fit1$loglik[2]+log(fit1$nevent)*sum(fit1$df), TRUE)
				}
				p <- which.min(myaic)
				ndf1 <- vector(length=nnl)
				for (s in 1: nnl) ndf1[s] <- out[p,s]
				aic <- myaic[p]
				ndf2 <- df2
				for (j in 1:nnl) {
					if (ndf1[j] == df2[j,1]) {
						a1 <- df2[j,1]
						a2 <- df2[j,2]
						a3 <- (df2[j,1]+df2[j,2])/2
						ndf2[j,] <- sort( c(a1, a2, a3) )
					}
					if (ndf1[j] == df2[j,3]) {
						a1 <- df2[j,3]
						a2 <- df2[j,2]
						a3 <- (df2[j,3]+df2[j,2])/2
						ndf2[j,] <- sort( c(a1, a2, a3) )
					}
					if (ndf1[j] == df2[j,2]) {
						p <- c()
						for ( t in 1:(3^nnl) ) {
							if ( out[t,j] != ndf1[j] & all(out[t,-j]-ndf1[-j] == 0) ) p <- c(p, t)
						}
						p1 <- which.min(myaic[p])
						if (p1 == 1) a2 <- df2[j,1]
						else a2 <- df2[j,3]
						a1 <- ndf1[j]
						a2 <- a2
						a3 <- (ndf1[j]+a2)/2
						a4 <- sort( c(a1, a2, a3) )
						ndf2[j,] <- a4
					}
				}
				df2<-ndf2
			}
		}
		ndf1 <- round(ndf1, 1)
		npred1 <- paste("pspline(", nl.predictors, ",df=", ndf1, ")", collapse="+")
		if ( !missing(other.predictors) ) pred1 <- paste(npred1, opred, sep="+")
		else pred1 <- npred1
		if (ctype == "TRUE") {
			covar <- as.formula( paste(" Surv(ntime,ntime2,nstatus)~ ", pred1) )
			try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
		} else {
			covar <- as.formula( paste(" Surv(ntime,nstatus)~ ", pred1) )
			try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
		}
		mydf <- fit1$df[1:nnl]
		aic <- -2*fit1$loglik[2]+2*sum(fit1$df)
		aicc <- -2*fit1$loglik[2]+2*( fit1$nevent*(sum(fit1$df)+1)/(fit1$nevent-sum(fit1$df)-2) )
		bic <- -2*fit1$loglik[2]+log(fit1$nevent)*sum(fit1$df)
	} else if (smoother == "ns") {
		if ( is.null(time2) ) mydf <- dfnsmult(time=time, status=status, nl.predictors=nl.predictors, other.predictors=other.predictors, method=method, mindf=mindf, maxdf=maxdf, ntimes=ntimes, data=data)
		else mydf <- dfnsmult(time=time, time2=time2, status=status, nl.predictors=nl.predictors, other.predictors=other.predictors, method=method, mindf=mindf, maxdf=maxdf, ntimes=ntimes, data=data)
		nnl1 <- seq(1:nnl)
		#ndf <- paste("df[", nnl1, "]", sep="")
		npred1 <- paste("ns(", nl.predictors, ",df=", mydf, ")", collapse="+")
		if ( !missing(other.predictors) ) pred1 <- paste(npred1, opred, sep="+")
		else pred1 <- npred1
		if (ctype == "TRUE") {
			p4 <- match(names(data), time2, nomatch=0)
			p5 <- which(p4 == 1)
			ntime2 <- data[,p5]
		}
		if (ctype == "TRUE") {
			covar <- as.formula( paste(" Surv(ntime,ntime2,nstatus)~ ", pred1) )
			fit1 <- coxph(covar, data=data, x=TRUE)
		} else {
			covar <- as.formula( paste(" Surv(ntime,nstatus)~ ", pred1) )
			fit1 <- coxph(covar, data=data, x=TRUE)
		}
		if ( !missing(other.predictors) ){
			nop <- length(other.predictors)
			mydf2 <- sum(mydf)+nop
			aic <- -2*fit1$loglik[2]+2*mydf2
			aicc <- -2*fit1$loglik[2]+2*( fit1$nevent*(mydf2+1)/(fit1$nevent-mydf2-2) )
			bic <- -2*fit1$loglik[2]+log(fit1$nevent)*mydf2
			#myaic <- -2*fit1$loglik[2]+2*sum(mydf)+2*nop
		} else {
			aic <- -2*fit1$loglik[2]+2*sum(mydf)
			aicc <- -2*fit1$loglik[2]+2*( fit1$nevent*(sum(mydf)+1)/(fit1$nevent-sum(mydf)-2) )
			bic <- -2*fit1$loglik[2]+log(fit1$nevent)*sum(mydf)
			#myaic <- -2*fit1$loglik[2]+2*sum(mydf)
		}
	}
	options(warn=0)
	return( list(df=mydf, AIC=aic, AICc=aicc, BIC=bic,  myfit=fit1, method=method, nl.predictors=nl.predictors) )
	#return( list(df=mydf) )
}
