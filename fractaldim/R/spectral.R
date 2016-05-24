do.fd.estimate1d.spectral <- function (x, y, use.lags, method.name, data.size,
							   method.name.main='', yname = '',
                               plot.loglog=FALSE, nlags = "auto", 
                               plot.allpoints=FALSE, legend.type='s', ...,
                               debuglevel=0) {
# This function is used by the periodogram and dctII
	nr.lags <- length(use.lags)
    rawD <- get.rawFD.from.regression(x[use.lags], y[use.lags])
  
  	if (debuglevel>=5) {
    	cat("\nLog-log plot output:\n")
    	summary(rawD)
  	}
  	if (plot.loglog) {
    	if (plot.allpoints & (length(x) > nr.lags)) {
      		rawD$x <- x
      		rawD$y <- y
    	}
    	plot.args <- list(...)
    	if (!is.element('main', names(plot.args))) plot.args[['main']] <- method.name.main
    	if (!is.element('ylab', names(plot.args))) plot.args[['ylab']] <- yname
    	if (!is.element('xlab', names(plot.args))) plot.args[['xlab']] <- paste("log(omega=[", 
    					 round(exp(rawD$x[1]),2),":",
                         round(exp(rawD$x[length(rawD$x)]),2),"])", sep="")
        filled <- if(plot.allpoints && (length(x) > nr.lags)) use.lags else NULL
    	do.call('fd.plot.loglog', c(list(rawD=rawD, filled=filled), plot.args))
  	}
  	DD <- list(	fd = if(!is.na(rawD$alpha)) (rawD$alpha+5)/2.0 else NA,
    			scale = if (!is.na(rawD$intercept)) exp(rawD$intercept) else NA    			)

  	if (plot.loglog) fd.plot.legend(rawD, DD, legend.type=legend.type)
  	
   	invisible (createFractalDim(dim=1, methods=method.name, FD=DD,
                              window.size=data.size, data.dim=1, loglog=rawD))
}


fd.estim.periodogram <- function (data, plot.loglog=FALSE, nlags = "auto", ...) {
  	rn<-length(data)
  	m <- (rn-1)/2
  	if (nlags == "auto") {
    	nr.lags <- trunc(min(m/2, rn^(2/3))) # suggestion of Chan, Hall, Poskitt
  	} else {
    	if (nlags != "all") {
      		nr.lags <- trunc(min(m/2, as.numeric(nlags)))
    	} else nr.lags <- trunc(m/2)
  	}
  	Jall<-ComputeSemipgramFFT(data)
  	npoints <- length(Jall)
  	nr.lags <- min(nr.lags, npoints)
  	omega <- seq(from=2*pi, to=2*npoints*pi, by=2*pi)
  	return(do.fd.estimate1d.spectral(x=log(omega[1:npoints]), y=log(Jall), use.lags=1:nr.lags,
  			method.name='periodogram', data.size=rn, method.name.main="Periodogram", 
  			yname=expression(log(J(omega))),
  			plot.loglog=plot.loglog, ...))
}


ComputeSemipgramFFT <- function(X)
  {
    N <- length(X)
    temp <- 2*(1:floor(floor(N/2)/2)) + 1
    return(((X[1] + X[N] + 2*Re(fft(c(0,X[2:(N-1)])))[temp])/(N-1))^2)
  }


fd.estim.dctII <- function (data, plot.loglog=FALSE, nlags = "auto", ...) { 
	rn<-length(data)
	nr.lags <- rn-1 
	if (nlags == "auto") { 
		nr.lags<- min(trunc(4*rn^(2/3)), nr.lags) # suggestion of Chan, Hall, Poskitt 
	} else { 
		if (nlags != "all") { 
			nlags <- as.numeric(nlags) 
			nr.lags <- min(nlags, nr.lags) 
		} 
	} 
	A <- ComputeDCT2FFT(data)[-1] 
	nr.lags <- min(nr.lags, length(A)) 
	logx<-log((1:length(A))*pi*(rn-1)/rn) 
	logAA<-log(A*A)
	return(do.fd.estimate1d.spectral(x=logx, y=logAA, use.lags=1:nr.lags, 
			method.name='dctII', data.size=rn, method.name.main="DCT-II", yname=expression(log(C2)), 
			plot.loglog=plot.loglog, ...)) 
} 

ComputeDCT2FFT <- function(X)
  {
    N <- length(X)
    return(Re(fft(c(X,rev(X)))[1:N] * 0.5 *
           exp(complex(N,rep(0,N),-pi*(0:(N-1))/(2*N))) *
           sqrt(c(1,rep(2,N-1))/N)))
  }

fd.estim.wavelet <- function (data,  plot.loglog=FALSE,
                           plot.allpoints=FALSE,
                           filter = "haar", 
                           J1=max(1,floor(log2(length(data))/3-1)), 
                           J0 = floor(log2(length(data))),
                           legend.type='s', ...,
                           debuglevel=0) {
  	require(wavelets)

  	result <- WaveVarFD(data, filter=filter, J1=J1, J0=J0, all.points=plot.allpoints)
  	rawD <- structure(list(alpha=result$beta, intercept=result$zeta,
                         	x=result$x[result$used], y=result$y[result$used], n=length(result$used), 
                         			lsq=NA),
                    class="FDloglog")
  	if (debuglevel>=5) {
    	cat("\nLog-log plot output:\n")
    	summary(rawD)
 	 }
  	if (plot.loglog) {
  		plot.args <- list(...)
    	if (!is.element('main', names(plot.args))) plot.args[['main']] <- "Wavelet"
    	if (!is.element('ylab', names(plot.args))) plot.args[['ylab']] <- "Y"
    	if (!is.element('xlab', names(plot.args))) plot.args[['xlab']] <- "log(tau)"
		if(plot.allpoints) {
    		rawD$x <- result$x
    		rawD$y <- result$y
    	}
    	filled <- if(plot.allpoints) result$used else NULL
    	do.call('fd.plot.loglog', c(list(rawD=rawD, filled=filled), plot.args))
  	}
	DD <- list(fd=result$fd, scale=exp(rawD$intercept))
	
  	if (plot.loglog) fd.plot.legend(rawD, DD, legend.type=legend.type)
  	
  	invisible (createFractalDim(dim=1, methods="wavelet", FD=DD,
                              window.size=length(data), data.dim=1, loglog=rawD))
}

WaveVarFD <- function(X, filter="haar", J1=max(1,floor(log2(length(X))/3-1)), J0 = floor(log2(length(X))), 
						all.points=FALSE)
  {
    N  <- length(X)
    J0.max <- floor(log2(N))
    J0.all <- if(all.points && J0 < J0.max) J0.max else min(J0.max, J0)
    Js <- 1:J0.all
    tau.all <- 2**(Js-1)
    log.tau.all <- log(tau.all)
    etao2 <- N/(2**(Js+1))  # used to be etao2 <- N/(2**Js)   
    N2 <- 2*N
    if (!is.element('numeric', class(X))) class(X) <- c('numeric', class(X)) # for modwt to work
    wt <- modwt(X, filter, J0.all, "reflection")
    wvar <- apply(sapply(wt@W, FUN=function(x) {x^2})[,Js] , MARGIN=2,
                  FUN=sum)/N2
    Y.all <- log(etao2*wvar) - digamma(etao2)
    w.all <- 1/trigamma(etao2)
    wxlog.tau.all <- w.all * log.tau.all
    Y <- Y.all[J1:J0]
    w <- w.all[J1:J0]
    log.tau <- log.tau.all[J1:J0]
    wxlog.tau <- w*log.tau.all[J1:J0]
    wxlog.tau2 <- w*((log.tau.all[J1:J0])^2)
    sum.w <- sum(w)
    sum.wxlog.tau <- sum(wxlog.tau)
    sum.wxlog.tau2 <- sum(wxlog.tau2)
    sum.wxY <- sum(w*Y)
    sum.wxlog.tauxY <- sum(wxlog.tau*Y)
    sum.wxlog.tau2 <- sum(w*(log.tau^2))
    bot <- sum.w*sum.wxlog.tau2 - sum.wxlog.tau^2
    ## intercept
    zeta <- (sum.wxlog.tau2*sum.wxY - sum.wxlog.tau*sum.wxlog.tauxY)/bot
    ## slope
    beta <- (sum.w*sum.wxlog.tauxY - sum.wxlog.tau*sum.wxY)/bot
    
    #beta <- (sum.wj*sum(wj.log.tauj*Y) - sum.wj.log.tauj*sum(Y))/
    #  (sum.wj*sum(wj.log.tauj^2) - sum.wj.log.tauj^2)
    #print(c("beta:", beta))
    #print(wj.log.tauj)
    #print(c("fd from beta1:", 2-beta1/2))

    return(list(zeta=zeta, beta=beta, fd=2 - beta/2, x=log.tau.all, y=Y.all, used=Js[J1:J0]))
  }
