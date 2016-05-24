dfnsuniv <- function(time, time2=NULL, status, nl.predictor, other.predictors, method, maxdf=NULL, data) {
	options(warn=-1)
	ctype <- "FALSE"
	if ( missing(data) ) stop("The argument data is missing")
	if ( missing(time) ) stop("The argument time is missing")
	if ( missing(status) ) stop("The argument status is missing")
	if ( missing(nl.predictor) ) stop("The argument 'nl.predictor' is missing")
	if ( !missing(time2) ) ctype <- "TRUE"
	if ( missing(method) ) method <- "AIC"
	if ( missing(maxdf) ) maxdf <- 15
	if ( method != "AIC" & method != "AICc" & method != "BIC") stop("The argument 'method' is not valid")
	p0 <- match(names(data), time, nomatch=0)
	p1 <- which(p0==1)
	ntime <- data[,p1]
	p2 <- match(names(data), status, nomatch=0)
	p3 <- which(p2==1)
	nstatus <- data[,p3]
	if (ctype == TRUE) {
		p2 <- match(names(data), time2, nomatch=0)
		p3 <- which(p2 == 1)
		ntime2 <- data[,p3]
	}
	myaic <- vector(length=maxdf)
	if ( !missing(other.predictors) ) opred <- paste(other.predictors, collapse="+")
	if (ctype == "TRUE") {
		p4 <- match(names(data), time2, nomatch=0)
		p5 <- which(p4 == 1)
		ntime2 <- data[,p5]
	}
	for ( k in 1:as.integer(maxdf) ) {
		npred1 <- paste("ns(", nl.predictor, ",df=", k, ")")
		if ( !missing(other.predictors) ) pred <- paste(npred1, opred, sep="+")
		else pred <- npred1
		#pred <- paste("ns(", predictor, ",df=", k, ")")
		if (ctype == "TRUE") {
			covar <- as.formula( paste(" Surv(ntime,ntime2,nstatus)~ ", pred) )
			try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
		} else {
			covar <- as.formula( paste(" Surv(ntime,nstatus)~ ", pred) )
			try(fit1 <- coxph(covar, data=data, x=TRUE), TRUE)
		}
		if (method == "AIC") try(myaic[k] <- -2*fit1$loglik[2]+2*k, TRUE)
		else if (method == "AICc") try(myaic[k] <- -2*fit1$loglik[2]+2*fit1$nevent*(k+1)/(fit1$nevent-k-2), TRUE)
		else if (method == "BIC") try(myaic[k] <- -2*fit1$loglik[2]+2*log(fit1$nevent)*k, TRUE)
	}
	df <- which.min(myaic)
	options(warn=0)
	return(df)
}
