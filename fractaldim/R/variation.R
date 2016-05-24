
do.fd.estimate1d.method <- function (data, fd.function, method.name, 
						  fd.fun.args = NULL, method.name.main='', yname = '', p = 2,
                          plot.loglog=FALSE, nlags = "auto", 
						  plot.allpoints=FALSE, legend.type='s', ...,
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
      		nr.lags <- min(length(x[x <= nlags]), n-1)
    	}
  	}  
  	#compute values of p-variation
  	g<-do.call(fd.function, c(list(data=data, lags=1:nr.lags, p=p), fd.fun.args))
  	x<-x[1:nr.lags]
  	logx<-log(x)
  	logg <- log(g)

    rawD <- get.rawFD.from.regression(logx, logg)

  	if (debuglevel>=5) {
    	cat("\nLog-log plot output:\n")
    	summary(rawD)
  	}
  	if (plot.loglog) {
  		grest <- NULL
    	if (plot.allpoints & (n-1) > nr.lags) {
      		grest <- rep(NA, n-1-nr.lags)
      		grest <- do.call(fd.function, c(list(data=data, lags=(1:(n-1-nr.lags)) + nr.lags, p=p),
      										 fd.fun.args))
      		logx <- log(1:(n-1))
      		logg <- c(log(g), log(grest))
      		rawD$x <- logx
      		rawD$y <- logg
    	}
    	plot.args <- list(...)
    	if (!is.element('main', names(plot.args))) plot.args[['main']] <- method.name.main
    	if (!is.element('ylab', names(plot.args))) plot.args[['ylab']] <- yname
    	if (!is.element('xlab', names(plot.args))) plot.args[['xlab']] <- paste("log(lags = [1:",
                         round((1:(n-1))[length(rawD$x)],2),"])", sep="")
		filled <- if(plot.allpoints && !is.null(grest)) 1:nr.lags else NULL
		do.call('fd.plot.loglog', c(list(rawD=rawD, filled=filled), plot.args))
  	}
  	DD <- list(	fd= if (!is.na(rawD$alpha)) 2-rawD$alpha/p else NA,
  				# normalized intercept
  				scale=if(!is.na(rawD$intercept)) exp(rawD$intercept)^(1/p) else NA
				)
  
  	if (plot.loglog) fd.plot.legend(rawD, DD, legend.type=legend.type)
  	
  	invisible (createFractalDim(dim=1, methods=method.name, FD=DD,
                              window.size=n, data.dim=1, loglog=rawD))
}

fd.estim.variation <- function (data, p.index=1, ...) {
	return(do.fd.estimate1d.method(data, 'fd.Variation', method.name="variation",
									method.name.main=paste("Variation - p =", p.index), 
									p=p.index, yname = "log(p-variation)", ...))
}

fd.estim.variogram <- function (data, ...) {
	p.index <- 2
	return(do.fd.estimate1d.method(data, 'fd.Variation', method.name="variogram",
									method.name.main="Variogram", 
									p=p.index, yname = "log(variogram)", ...))
}

fd.estim.madogram <- function (data, ...) {
	p.index <- 1
	return(do.fd.estimate1d.method(data, 'fd.Variation', method.name="madogram",
									method.name.main="Madogram", 
									p=p.index, yname = "log(madogram)", ...))
}

fd.estim.rodogram <- function (data, ...) {
	p.index <- 0.5
	return(do.fd.estimate1d.method(data, 'fd.Variation', method.name="rodogram",
									method.name.main="Rodogram", 
									p=p.index, yname = "log(rodogram)", ...))
}


fd.Variation <- function(data, lags=1:2, p=1, ...) {
  	n <- length(data)
  	m <- length(lags)
  	g<-rep(NA, m)
  	for (j in 1:m) {
  		lag <- lags[j]
    	g[j] <- sum(abs(data[(lag+1):n]-data[1:(n-lag)])^p)/(2*(n-lag))
  	}
  	return(g)
}


fd.estim.incr1 <- function (data, p.index=2, ...) {
	return(do.fd.estimate1d.method(data, 'fd.Increment1', method.name="incr1",
									method.name.main=paste("Increment1 - p =", p.index), 
									p=p.index, yname = "log(p-variation)", ...))
}

fd.Increment1 <- function(data, lags=1:2, p=2, ...) {
	n <- length(data)
  	m <- length(lags)
  	g<-rep(NA, m)
  	for (j in 1:m) {
  		lag <- lags[j]
  		if (lag > n/2-1) break
    	g[j] <- sum(abs(data[1:(n-2*lag)]-2*data[(lag+1):(n-lag)]+data[(2*lag+1):n])^p)/(n-2*lag)
  	}
  	return(g)
}

fd.estim.genton <- function (data, ...) {
	require(pcaPP)
	return(do.fd.estimate1d.method(data, 'fd.Genton', method.name='genton',
									method.name.main="Robust Genton", 
									p=2, yname = "log(Y)", ...))
}

fd.Genton <- function(data, lags=1:2, ...) {
	n <- length(data)
  	m <- length(lags)
  	g<-rep(NA, m)
  	for (j in 1:m) {
  		lag <- lags[j]
  		Vh <- data[(lag+1):n]-data[1:(n-lag)]
  		g[j] <- (qn(Vh))^2
  	}
	return(g)
}