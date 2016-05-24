smoothHR <- function(data, time=NULL, time2=NULL, status=NULL, formula=NULL, coxfit, status.event=NULL) {
	ctype <- "FALSE"
	modelfit <- "TRUE"
	mydata2 <- deparse( substitute(data) )
	if ( missing(coxfit) ) modelfit <- "FALSE"
	if ( !missing(coxfit) ) if ( !inherits(coxfit, "coxph") ) stop("Argument coxfit must be of class coxph")
	if ( !is.data.frame(data) ) stop("data must be of class data.frame")
	if (modelfit == "TRUE") {
		if ( !inherits(coxfit, "coxph") ) stop("Object coxfit must be of class coxph")
		if ( is.null(coxfit$x) ) stop("The argumment x in the coxph object is missing")
			fit <- coxfit
	}
	#if (!missing(formula) & !missing(coxfit)) stop("....only one is requested")
	if ( missing(data) ) stop("The argumment data is missing")
	if (missing(time) & modelfit == "FALSE") stop("The argumment time is missing")
	if (missing(status) & modelfit == "FALSE") stop("The argumment status is missing")
	if (missing(formula) & modelfit == "FALSE") stop("The argumment formula is missing")
	if (!missing(time2) & modelfit == "FALSE") ctype <- "TRUE"
	mydata <- data
	if (modelfit == "FALSE") {
		p0 <- match(names(data),time, nomatch=0)
		p1 <- which(p0 == 1)
		ntime <- data[,p1]
		p2 <- match(names(data), status, nomatch=0)
		p3 <- which(p2 == 1)
		nstatus <- data[,p3]
		if (ctype == "TRUE") {
			p4 <- match(names(data), time2, nomatch=0)
			p5 <- which(p4 == 1)
			ntime2 <- data[,p5]
		}
		if ( !missing(time) ) time <- ntime
        if ( !missing(time2) ) time2 <- ntime2
        if ( !missing(status) ) status <- nstatus
        #if (!missing(status) & (nlevels(as.factor(nstatus))!=2 | (levels(as.factor(nstatus))[1]!="0" | levels(as.factor(nstatus))[2]!="1")))  stop("status must be of class factor with levels '0' and '1' with '1' for the event")
        if ( !missing(time) & !is.numeric(ntime) & !is.integer(ntime) ) stop("time must be of class numeric or integer")
        if ( !missing(time2) ) if ( !is.numeric(ntime2) & !is.integer(ntime2) ) stop("time2 must be of class numeric or integer")
        if ( !missing(status) & missing(status.event) ) status.event <- max(nstatus)
		fmla <- attr(terms(formula), "term.labels")
		ncov <- length(fmla)
		colvar <- rep(0, ncov)
		for (k in 1:ncov) {
			if ( fmla[k] %in% names(data) ) {
				colvar[k] <- which(names(data) == fmla[k])
            } else {
                for ( j in 1:ncol(data) ) {
					if ( any( grep(names(data)[j], fmla[k]) ) ) colvar[k] <- j
				}
			}
		}
		if ( any(colvar == 0) ) stop("'formula' must contain the right variables")
		if (ctype == "TRUE") {
			covar <- as.formula( paste( " Surv(ntime,ntime2,nstatus)~ ", paste(fmla, collapse = "+") ) )
            covar2 <- as.formula( paste( " Surv(ntime,ntime2,nstatus==", status.event, ")~ ", paste(fmla, collapse = "+") ) )
			fit <- coxph(covar2, data = data, x=T)
		} else {
			covar <- as.formula( paste( " Surv(ntime,nstatus)~ ", paste(fmla, collapse = "+") ) )
			covar2 <- as.formula( paste(" Surv(ntime,nstatus==", status.event, ")~ ", paste(fmla, collapse = "+") ) )
			fit <- coxph(covar2, data=data, x=TRUE)
		}
	}
	phtest <- cox.zph(fit)
	dimt <- dim(phtest$table)
	p_value <- phtest$table[dimt[1],dimt[2]]
	if (p_value <= 0.05) cat("Warning: a global p-value of ", p_value, "was obtained when testing the proportional hazards assumption", "\n")
	a1 <- c()
	if (modelfit == "TRUE") {
		for(k in 1:dim(mydata)[2]) {
			if ( any( grep(names(mydata)[k], fit$call) ) ) a1 <- c(a1,k)
		}
		mydata <- mydata[,a1]
	}
	if (modelfit == "FALSE") {
		if (ctype == "TRUE") a1 <- c(p1, p5, p3)
		else a1 <- c(p1, p3)
		for (k in 1:dim(mydata)[2]) {
			if ( any( grep(names(mydata)[k], formula) ) ) a1<-c(a1,k)
		}
		mydata <- mydata[,a1]
	}
	if (modelfit == "TRUE") if (as.name( toString(fit$call[[3]]) ) != mydata2) cat("Warning: check if argument 'data' is the same used in coxfit", "\n")
	mydata <- na.omit(mydata)
	nv <- fit$nevent
	if (ctype == "TRUE") {
		myd <- as.data.frame( table(mydata[,3]) )
		w.status.event <- which(myd[,1] == status.event)
		if (nv != table( as.factor(mydata[,3]) )[[w.status.event]]) cat("Warning: check if argument 'data' is the same used in coxfit", "\n")
	}
	object <- list(dataset=mydata, coxfit=fit, phtest=phtest)
	class(object) <- "HR"
	return(object)
}
