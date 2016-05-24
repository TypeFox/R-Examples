## Copyright (C) 2006-2008  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

#' @export
llar <- function(x, m, d=1, steps=d, series, eps.min=sd(x)/2, eps.max=diff(range(x)), neps=30, trace=0) {
	if(missing(series))
		series <- deparse(substitute(x))
  if(NCOL(x)>1)
    stop("Multivariate time series not allowed")
  if(any(is.na(x)))
    stop("NA's not allowed")
	epsSeq <- exp(seq(log(eps.min), log(eps.max), length=neps))
	epsSeq <- epsSeq/diff(range(x))
	x1 <- (x - min(x))/(diff(range(x)))
	err <- rep(-1, length(epsSeq))
        nok <- rep(0, length(epsSeq))
        avfound <- rep(-1, length(epsSeq))
	res <- .C("llar", series=as.double(x1), length=as.integer(length(x1)),
                  m=as.integer(m), d=as.integer(d), steps=as.integer(steps), tol=as.double(1e-7),
                  epsSeq = as.double(epsSeq), NEPS=length(epsSeq),
                  trace=as.integer(trace), err=as.double(err), nok=as.integer(nok),
                  avfound=as.double(avfound), PACKAGE="tsDyn")
        err <- res$err
        err[err==-1] <- NA
        nok <- res$nok
        avfound <- res$avfound
        avfound[avfound==-1] <- NA
	res <- structure(list(RMSE=err, eps=epsSeq*diff(range(x)),
		frac = nok/(length(x)-(m-1)*d), avfound = avfound), 
		series=series, call=match.call(), class="llar")
	return( res )
}

#' @S3method print llar
print.llar <- function(x, ...) {
	cat("\nCall: ")
	print(attr(x,"call"))
	reps <- range(x$eps)
	cat("\nNeighborhood size ranging from ", reps[1]," to ", reps[2],"\n")
	reps <- range(x$RMSE)
	cat("Relative error ranging from ", reps[1]," to ", reps[2],"\n\n")
}

#' @S3method plot llar
plot.llar <- function(x, ...)
	plot(x$eps, x$RMSE, xlab=expression(epsilon), ylab="relative error", log="x", type="l",
		ylim=c(0,max(1,max(x$RMSE))),
		main=paste("local linear fit of",attr(x,"series")))

#' @S3method as.data.frame llar
as.data.frame.llar <- function(x, row.names, optional, ...)
	data.frame(unclass(x))

#One steps ahead forecast from local linear model applied to the observed time series.
llar.step <- function(x, m, d=1, steps=d, series, eps=stop("you must specify a window value"),
                         onvoid=c("fail","enlarge"), r = 20, trace=1) {
	if(missing(series))
		series <- deparse(substitute(x))
  if(NCOL(x)>1)
    stop("Multivariate time series not allowed")
  if(any(is.na(x)))
    stop("NA's not allowed")  
  onvoid <- match.arg(onvoid)
  r <- 1+r/100
  lags <- c(steps, 0:(m-1)*(-d))
  xxFull <- embedd(x, lags=lags)
  yy <- xxFull[,1]
  xx <- xxFull[,seq_len(m)+1]
  xxL <- xx[nrow(xx),]
  dif <- sweep(xx, 2, xxL, "-")
  DL <- apply(dif,1,function(x) sum(x^2))
  DL <- sqrt(DL)#euclidean norm
  SEL <- which(DL<=eps)
  SEL <- SEL[abs(SEL-nrow(xx))>=steps]#exclude theiler window
  pfound <- length(SEL)
  if(pfound<(2*m+1)){
    if(onvoid=="fail")
      stop("not enough neighbours found")
    else {
      while(pfound<(2*m+1)) {
        eps <- eps * r
        SEL <- which(DL<=eps)
        SEL <- SEL[abs(SEL-nrow(xx))>=steps]#exclude theiler window
        pfound <- length(SEL)
      }
      if(trace)
        cat("Isolated point. Using eps =",format(eps,digits=4),"\n")
    }
  }
  xxSEL <- cbind(1,xx[SEL,])
  yySEL <- yy[SEL]
  lmod <- lm.fit(xxSEL, yySEL)
  return(crossprod(c(1,xxFull[nrow(xxFull), seq_len(m)]),lmod$coef))
}

#' @export
llar.predict <- function(x, m, d=1, steps=d, series, n.ahead=1,
                         eps=stop("you must specify a window value"),
                         onvoid=c("fail","enlarge"), r = 20, trace=1) {
	if(missing(series))
		series <- deparse(substitute(x))
  if(NCOL(x)>1)
    stop("Multivariate time series not allowed")
  if(any(is.na(x)))
    stop("NA's not allowed")  
  n <- length(x)
  for(i in 1:n.ahead) {
    x <- c(x, llar.step(x, m=m, d=d, steps=steps, series=series, eps=eps,onvoid=onvoid,r=r,trace=trace))
  }
  return(x[(n+1):(n+n.ahead)])
}

#' @export
llar.fitted <- function(x, m, d=1, steps=d, series, eps, trace=0) {
	if(missing(series))
		series <- deparse(substitute(x))
  if(NCOL(x)>1)
    stop("Multivariate time series not allowed")
  if(any(is.na(x)))
    stop("NA's not allowed")
  x1 <- (x- min(x))/(diff(range(x)))
  blength <- length(x)-(m-1)*d -steps
  fitted <- rep(0, blength)
  nOK <- numeric(blength)
  res <- .C("fittedllar", series=as.double(x1), length=as.integer(length(x1)), 
            m=as.integer(m), d=as.integer(d), steps=as.integer(steps), tol=as.double(1e-7),
            eps = as.double(eps), trace=as.integer(trace),
            fitted=as.double(fitted), nOK=as.integer(nOK), PACKAGE="tsDyn")
  mask <- (res$nOK>(2*(m+1)) )
  fitted[mask] <- res$fitted
  fitted[!mask] <- NA
  fitted <- fitted * diff(range(x)) + min(x)
  return( fitted )
}
