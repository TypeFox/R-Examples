fd.estim.boxcount <- function (data, plot.loglog=FALSE, nlags = "auto", 
							shift.up=TRUE, plot.allpoints=FALSE, 
							legend.type='s', ...,
                            debuglevel=0) {
  	n <- length(data)
	m <- trunc(log2(n-1))+1
  	neff<-trunc((n-1)/(2^(m-1)))*(2^(m-1))+1 # n-1 dyadic

  	nr.lags <- m  
	if (nlags == "auto") {
    	nr.lags<-max(2,nr.lags-2) # exclude the two largest boxes (Liebovitch and Toth)
  	} else {
    	if (nlags != "all") {
      		nlags <- as.numeric(nlags)
      		nr.lags <- min(nlags, nr.lags)
    	}
  	}
	data.range <- range(data[1:neff])
	N <- rep(NA, m)
  	sizes <- 2^(0:(m-1))
  	for (ll in nr.lags:1) { # iterate over box sizes
    	N[ll] <- ComputePlotPoint.boxcount(data, neff, sizes[ll], 
    							data.range=data.range, shift.up=shift.up)
    	if ((nlags == "auto") && (N[ll] > (neff/5)) && (ll < (nr.lags-1))) { # Liebovitch and Toth modification 
      		N[ll]<-NA
      		break
    	}
  	}
  	nin <- !is.na(N)
  	logx<-log(sizes[nin])
  	logN<-log(N[nin])

    rawD <- get.rawFD.from.regression(logx, logN)
  
 	if (debuglevel>=5) {
    	cat("\nLog-log plot output:\n")
    	summary(rawD)
  	}
  	if (plot.loglog) {
    	if (plot.allpoints) {
      		isnan <- is.na(N)
      		for (ll in (1:m)[isnan]) {
				N[ll] <- ComputePlotPoint.boxcount (data, neff, sizes[ll], 
									data.range=data.range, shift.up=shift.up)
      		}
      		logx <- log(sizes)
      		logN <- log(N)
      		rawD$x <- logx
      		rawD$y <- logN
    	}
    	l <- 0:(length(rawD$x)-1)
    	plot.args <- list(...)
    	if (!is.element('main', names(plot.args))) plot.args[['main']] <- "Box-Count"
    	if (!is.element('ylab', names(plot.args))) plot.args[['ylab']] <- "log(N(l))"
    	if (!is.element('xlab', names(plot.args))) plot.args[['xlab']] <- paste(
    				"log(2^l), l = [", l[1],":",l[length(l)],"]", sep="")

		filled <- if(plot.allpoints && sum(isnan) < length(N)) nin else NULL
    	do.call('fd.plot.loglog', c(list(rawD=rawD, filled=filled), plot.args))
  	}
	DD <- list(	fd=if (!is.na(rawD$alpha)) - rawD$alpha else NA,
			  	scale=if (!is.na(rawD$intercept)) exp(rawD$intercept) else NA
			  )
  
  	if (plot.loglog) fd.plot.legend(rawD, DD, legend.type=legend.type)
	
  	invisible (createFractalDim(dim=1, methods=c("boxcount"), FD=DD,
                              window.size=neff, data.dim=1, loglog=rawD))
}

ComputePlotPoint.boxcount <- function(data, n, boxsize, data.range, shift.up=TRUE) {
	# 'shift.up' causes shifting up blocks vertically
	N <- 0
  	width <- boxsize/n
  	height <- diff(data.range)*boxsize/(n-1)
  	minblock <- data.range[1]
  	qh <- as.integer((n-1)/boxsize)
  	qv <- trunc(diff(data.range)/height)
  	for (hor in 1:qh) { # iterate over x-axis
    	block <- data[((hor-1)*boxsize+1):(hor*boxsize+1)]
    	nblock <- length(block)
    	maxblock <-max(block)
    	if (shift.up) minblock <- min(block)
    	for (ver in 1:qv) { # iterate over the y-axis
    		box <- c(minblock+(ver-1)*height, minblock+ver*height)
      		if (maxblock < box[1]) break
      		if (minblock > box[2]) break
        	N <- if(IsInBox(block,box,nblock)) N+1 else N + CrossingBox(block,box,nblock)
      	}
    }
  	return(N)
}

IsInBox <- function(data,box,ld) {
	# Is any of the data points in the box?
	# box is a 2-elements vector of (low, up)
	if(any((data >= box[1]) & (data <= box[2]))) return(TRUE)
	return (FALSE)
}

CrossingBox <- function(data,box,ld) {
	# is the given box crossed by lines connecting points outside of the box? 
	# box is a 2-elements vector of (low, up)
	dgt.up <- (1:ld)[data > box[2]]
	dlt.low <- (1:ld)[data < box[1]]
	# find neighbours for which the two conditions above hold
	if(any(is.element(dlt.low+1,dgt.up))) return(TRUE)
	if(any(is.element(dlt.low,dgt.up+1))) return(TRUE)
  	return(FALSE)
}
