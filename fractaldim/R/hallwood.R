fd.estim.hallwood <- function (data, plot.loglog=FALSE, nlags = "auto", 
								plot.allpoints=FALSE, 
                                legend.type='s', ...,
                                debuglevel=0) {
	n<-length(data)
  	x<-1:n

  	if (nlags == "auto") {
    	nr.lags <- 2
  	} else {
    	if (nlags == "all") {
      		nr.lags <- n-1
    	} else {
      		nlags <- as.numeric(nlags)
      		nr.lags <- min(length(x[x<= nlags]), n-1)
    	}
  	}
  	a <- fd.HallWood(data, k=nr.lags, FALSE)
  	x<-x[1:nr.lags]
  	logx<-log(x)
  	loga <- log(a)

    rawD <- get.rawFD.from.regression(logx, loga)
  	
  	if (debuglevel>=5) {
    	cat("\nLog-log plot output:\n")
    	summary(rawD)
  	}
  	if (plot.loglog) {
  		arest <- NULL
    	if (plot.allpoints & (n-1) > nr.lags) {
      		arest <- rep(NA, n-1-nr.lags)
      		for (ll in 1:(n-1-nr.lags)) {
        		arest[ll] <- ComputePlotPoint.hallwood(data, n, ll+nr.lags, FALSE)
      		}
      		logx <- log(1:(n-1))
      		loga <- c(log(a), log(arest))
      		rawD$x <- logx
      		rawD$y <- loga
    	}
    	plot.args <- list(...)
    	if (!is.element('main', names(plot.args))) plot.args[['main']] <- "Hall-Wood"
    	if (!is.element('ylab', names(plot.args))) plot.args[['ylab']] <- "log(A)"
    	if (!is.element('xlab', names(plot.args))) plot.args[['xlab']] <- paste("log(lags = [1:",
                         round((1:(n-1))[length(rawD$x)],2),"])", sep="")
        filled <- if(plot.allpoints && !is.null(arest)) 1:nr.lags else NULL
    	do.call('fd.plot.loglog', c(list(rawD=rawD, filled=filled), plot.args))
  	}
  	DD <- list(	fd = if (!is.na(rawD$alpha)) 2  - rawD$alpha else NA,
    			scale = if (!is.na(rawD$intercept)) exp(rawD$intercept) else NA)
				
  	if (plot.loglog) fd.plot.legend(rawD, DD, legend.type=legend.type)
  	
  	invisible (createFractalDim(dim=1, methods=c("hallwood"), FD=DD,
                              window.size=n, data.dim=1, loglog=rawD))
}

fd.HallWood <- function (data, k=2, neff=TRUE) {
  	n <- length(data)
  	a<-rep(NA, k)
  	for (l in 1:k) {
    	a[l] <- ComputePlotPoint.hallwood(data, n, l, neff)
  	}
  	return (a)
}

ComputePlotPoint.hallwood <- function (data, n, lag, neff=TRUE) {
  	q<-trunc((n-1)/(lag))
  	epsilon <- lag/n
  	mya <- 0
  	for (i in 1:q) {
    	block <- data[seq(from=(i-1)*lag+1, to=i*lag+1, by=lag)]
    	mya <- mya + (max(block) - min(block))
  	}
  	a <- epsilon*mya
  	return (ifelse(neff, a*(n-1)/(lag*q),a))
}
