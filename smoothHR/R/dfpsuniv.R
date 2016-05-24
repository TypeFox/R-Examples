dfpsuniv <- function(time, time2=NULL, status, nl.predictors, other.predictors, smoother, method, mindf=NULL, maxdf=NULL, ntimes=NULL, data) {
	options(warn=-1)
	ctype <- "FALSE"
	if ( missing(data) ) stop("The argument data is missing")
	if ( missing(time) ) stop("The argument time is missing")
	if ( missing(status) ) stop("The argument status is missing")
	if ( missing(nl.predictors) ) stop("The argument 'nl.predictors' is missing")
	if ( !missing(time2) ) ctype <- "TRUE"
	if ( missing(smoother) ) smoother <- "ns"
	if (smoother != "ns" & smoother != "pspline") stop("argument 'smoother' must be 'ns' or 'pspline'")
	if ( missing(method) ) method <- "AIC"
	if ( missing(mindf) ) mindf <- 1.5
	if ( missing(ntimes) ) ntimes <- 5
	if ( method != "AIC" & method != "AICc" & method != "BIC") stop("The argument 'method' is not valid")
	if (ntimes < 2 | ntimes > 10) stop("The argument 'ntimes' must be between 2 and 10")
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
	p2 <- match(names(data), nl.predictors[1], nomatch=0)
	if (sum(p2)==0) stop("Check variables in argument 'nl.predictors'")
	if ( !missing(other.predictors) ) {
		nop <- length(other.predictors)
		for (i in 1:nop) {
			p2 <- match(names(data), other.predictors[i], nomatch=0)
			if (sum(p2) == 0) stop("Check variables in argument 'other.predictors'")
		}
	}
	if (dim(table(nl.predictors)) < 1) stop("Check variables in argument nl.predictors")
	if ( !missing(other.predictors) ) {
		nop <- length(other.predictors)
		if (dim( table(other.predictors) ) < nop) stop("Check variables in argument other.predictors")
		all.predictors <- c(other.predictors, nl.predictors)
		if (dim( table(all.predictors) ) < nop+1) stop("Check variables in argument's nl.predictors and other.predictors")
	}
	if (method == "AIC") {
		if ( !missing(other.predictors) ) opred <- paste(other.predictors, collapse="+")
		npred0 <- paste("pspline(", nl.predictors, ", df=0)", collapse="+")
		if ( !missing(other.predictors) ) pred <- paste(npred0, opred, sep="+")
		else pred <- npred0
	} else if (method == "AICc") {
		if ( !missing(other.predictors) ) opred <- paste(other.predictors, collapse="+")
		npred0 <- paste("pspline(", nl.predictors, ", df=0, caic=TRUE)", collapse="+")
		if ( !missing(other.predictors) ) pred <- paste(npred0, opred, sep="+")
		else pred <- npred0
	} else if (method=="BIC") {
		if ( !missing(other.predictors) ) opred <- paste(other.predictors, collapse="+")
		npred0 <- paste("pspline(", nl.predictors, ", df=0)", collapse="+")
		if ( !missing(other.predictors) ) pred <- paste(npred0, opred, sep="+")
		else pred <- npred0
	}
	#ndf <- paste("df[",1,"]", sep="")
	if (ctype == "TRUE") {
		p4 <- match(names(data), time2, nomatch=0)
		p5 <- which(p4 == 1)
		ntime2 <- data[,p5]
		covar <- as.formula( paste(" Surv(ntime,ntime2,nstatus)~ ", pred) )
		fit <- coxph(covar, data=data, x=TRUE)
	} else {
		covar <- as.formula( paste(" Surv(ntime,nstatus)~ ", pred) )
		fit <- coxph(covar, data=data, x=TRUE)
	}
	fit1<-fit
	ndf1 <- fit$df[1]
	if ( !missing(maxdf) ) ndf1 <- maxdf
	if (method == "BIC") {
		df0 <- 1.5
		if ( !missing(mindf) ) df0 <- mindf
		df1 <- ndf1
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
				try(myaic[k] <- -2*fit1$loglik[2]+log(fit1$nevent)*sum(fit1$df), TRUE)
			}
			p <- which.min(myaic)
			ndf1 <- vector(length=1)
			ndf1 <- out[p,1]
			aic <- myaic[p]
			ndf2 <- df2
			if (ndf1 == df2[1,1]) {
				a1 <- df2[1,1]
				a2 <- df2[1,2]
				a3 <- (df2[1,1]+df2[1,2])/2
				ndf2[1,] <- sort( c(a1, a2, a3) )
			}
			if (ndf1 == df2[1,3]) {
				a1 <- df2[1,3]
				a2 <- df2[1,2]
				a3 <- (df2[1,3]+df2[1,2])/2
				ndf2[1,] <- sort( c(a1, a2, a3) )
			}
			if (ndf1 == df2[1,2]) {
				p <- c()
				for ( t in 1:3 ) {
					if ( out[t,1] != ndf1 ) p <- c(p, t)
				}
				p1 <- which.min(myaic[p])
				if (p1 == 1) a2 <- df2[1,1]
				else a2 <- df2[1,3]
				a1 <- ndf1
				a2 <- a2
				a3 <- (ndf1+a2)/2
				a4 <- sort( c(a1, a2, a3) )
				ndf2[1,] <- a4
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
	mydf <- fit1$df[1]
	options(warn=0)
	return(mydf)
}
