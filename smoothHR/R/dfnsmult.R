dfnsmult <- function(time, time2=NULL, status, nl.predictors, other.predictors, method, mindf=NULL, maxdf=NULL, ntimes=NULL, data) {
	options(warn=-1)
	ctype <- "FALSE"
	if ( missing(data) ) stop("The argument data is missing")
	if ( missing(time) ) stop("The argument time is missing")
	if ( missing(status) ) stop("The argument status is missing")
	if ( missing(nl.predictors) ) stop("The argument 'nl.predictors' is missing")
	if ( !missing(time2) ) ctype <- "TRUE"
	if ( missing(method) ) method <- "AIC"
	if ( missing(ntimes) ) ntimes <- 5
	if ( method != "AIC" & method != "AICc" & method != "BIC") stop("The argument 'method' is not valid")
	if (ntimes < 2 | ntimes > 10) stop("The argument 'ntimes' must be between 2 and 10")
	p0 <- match(names(data), time, nomatch=0)
	p1 <- which(p0 == 1)
	ntime <- data[,p1]
	p2 <- match(names(data), status, nomatch=0)
	p3 <- which(p2 == 1)
	nstatus <- data[,p3]
	if (ctype == TRUE) {
		p2 <- match(names(data), time2, nomatch=0)
		p3 <- which(p2 == 1)
		ntime2 <- data[,p3]
	}
	if ( !missing(other.predictors) ) opred <- paste(other.predictors, collapse="+")
	#if (!missing(other.predictors)) opred<-paste(other.predictors, collapse = "+")
	#npred0<-paste("ns(",nl.predictors,",df=0,caic=TRUE)",collapse="+")
	#if (!missing(other.predictors)) pred<-paste(npred0,opred,sep="+")
	#else pred<-npred0
	nnl <- length(nl.predictors)
	if ( missing(mindf) ) mindf <- rep(1,nnl)
	if ( missing(maxdf) ) maxdf <- rep(15,nnl)
	if ( nnl == 1 & is.null(time2) ) ndf1 <- dfnsuniv(time=time, status=status, nl.predictor=nl.predictors, other.predictors=other.predictors, method=method, maxdf=maxdf, data=data)
	if ( nnl == 1 & !is.null(time2) ) ndf1 <- dfnsuniv(time=time, time2=time2, status=status, nl.predictor=nl.predictors, other.predictors=other.predictors, method=method, maxdf=maxdf, data=data)
	if (nnl > 1) {
		if (nnl > 10) stop("The maximum number of nonlinear predictors is 10")
		if (dim( table(nl.predictors) ) < nnl) stop("Check variables in argument nl.predictors")
		if ( !missing(other.predictors) ) {
			nop <- length(other.predictors)
			if (dim( table(other.predictors) ) < nop) stop("Check variables in argument other.predictors")
			all.predictors <- c(other.predictors, nl.predictors)
			if (dim( table(all.predictors) ) < nop+nnl) stop("Check variables in argument's nl.predictors and other.predictors")
		}
		if (nnl == 5) cat("This may take a few seconds...\n")
		nnl1 <- seq(1:nnl)
		ndf <- paste("df[", nnl1, "]", sep="")
		df0 <- c( rep(1, nnl) )
		df1 <- c( rep(15, nnl) )
		if ( !missing(mindf) ) df0 <- mindf
		if ( !missing(maxdf) ) df1 <- maxdf
		df2 <- cbind(df0, round( (df0+df1)/2, 0 ), df1)
		for (i in 1:ntimes) {                                                                                    
			aaa <- lapply( apply( df2, 1, function(z) {list( c(z[1],z[2],z[3]) );} ), function(y) {unlist(y);} )
			out <- do.call(expand.grid, aaa)
			myaic <- rep(100000, dim(out)[1])
			for (k in 1:dim(out)[1]) {
				df <- round(out[k,], 0)
				npred1 <- paste("ns(", nl.predictors, ",df=", df, ")", collapse="+")
				if ( !missing(other.predictors) ) pred1 <- paste(npred1, opred, sep="+")
				else pred1 <- npred1
				if (ctype == "TRUE") {
					covar <- as.formula( paste(" Surv(ntime,ntime2,nstatus)~ ", pred1) )
					try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
				} else {
					covar <- as.formula( paste(" Surv(ntime,nstatus)~ ", pred1) )
					try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
				}
				#try(myaic[k] <- -2*fit1$loglik[2]+2*sum(out[k,]), TRUE)
				if (method == "AIC") try(myaic[k] <- -2*fit1$loglik[2]+2*sum(out[k,]), TRUE)
				else if (method == "AICc") try(myaic[k] <- -2*fit1$loglik[2]+2*fit1$nevent*(sum(out[k,])+1)/(fit1$nevent-sum(out[k,])-2), TRUE)
				else if (method == "BIC") try(myaic[k] <- -2*fit1$loglik[2]+2*log(fit1$nevent)*sum(out[k,]), TRUE)
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
			df2 <- ndf2
		}
    }
	options(warn=0)
	return( round(ndf1) )
}
