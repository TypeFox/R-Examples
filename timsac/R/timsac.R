.packageName <- "timsac"

mfilter <-
function (x, filter, method = c("convolution", "recursive"), init = NULL) 
{
    method <- match.arg(method)
    x <- as.ts(x)
    xtsp <- tsp(x)
    x <- as.matrix(x)
    n <- nrow(x)
    nser <- ncol(x)
    p <- dim(filter)[3]
    nd <- n+p
    if (any(is.na(filter))) 
        stop("missing values in 'filter'")
    y <- matrix(0, nd, nser)
    if (method == "convolution") {
        if (p > n) 
            stop("'filter' is longer than time series")
        if (missing(init)) {
            init <- matrix(0, p, nser)
        }
        else {
            ni <- NROW(init)
            if (ni != p) 
                stop("length of 'init' must equal length of 'filter'")
            if (NCOL(init) != 1 && NCOL(init) != nser) 
                stop(gettextf("'init'; must have 1 or %d cols", 
                  nser), domain = NA)
            if (!is.matrix(init)) 
                init <- matrix(init, p, nser)
        }
        xx <- matrix(NA,nd,nser)
        for (i in 1:p) xx[i, ] <- init[i, ]
        tar <- array(0,dim=c(p,nser,nser))
        for (i in 1:p) tar[i, ,] <- t(filter[, ,i])
        i <- p + 1
        while (i <= nd) {
          xx[i,] = x[i-p,]
          y[i, ] <- x[i-p, ]
          for (j in 1:p) y[i, ] <- y[i, ] - xx[i - j, ] %*% tar[j, ,]
             i <- i + 1
          }
    }
    else {
        if (missing(init)) {
            init <- matrix(0, p, nser)
        }
        else {
            ni <- NROW(init)
            if (ni != p) 
                stop("length of 'init' must equal length of 'filter'")
            if (NCOL(init) != 1 && NCOL(init) != nser) 
                stop(gettextf("'init'; must have 1 or %d cols", 
                  nser), domain = NA)
            if (!is.matrix(init)) 
                init <- matrix(init, p, nser)
        }
        for (i in 1:p) y[i, ] <- init[i, ]
        tar <- array(0,dim=c(p,nser,nser))
        for (i in 1:p) tar[i, ,] <- t(filter[, ,i])
        i <- p + 1
        while (i <= nd) {
          y[i, ] <- x[i-p, ]
          for (j in 1:p) y[i, ] <- y[i, ] + 
                     y[i - j, ] %*% tar[j, ,]
             i <- i + 1
          }
     }
    y <- y[(p+1):nd, ]
    y <- drop(y)
    tsp(y) <- xtsp
    class(y) <- "mts"
    return(y)
}

#####   TIMSAC72   #####

autcor <- function (y, lag=NULL, plot=TRUE, lag_axis=TRUE)
{
    n <- length(y)
    if( is.null(lag) )  lag <- as.integer(2*sqrt(n))    # maximum lag
    lag1 <- lag+1

    z <- .Call("autcor",
             as.double(y),
             as.integer(n),
             as.integer(lag1))

    if( plot == TRUE ) {
      plot((0:lag), z[[1L]], type="h", ylab="Autocorrelation", xlab="Lag")
      if( lag_axis == TRUE ) abline(h=0, lty=1)
    }

    autcor.out <- list( acov=z[[1L]], acor=z[[2L]], mean=z[[3L]] )
    return( autcor.out )
}


fftcor <- function (y, lag=NULL, isw=4, plot=TRUE, lag_axis=TRUE)
{
    ld <- nrow(y)                # length of data
    d <- ncol(y)                 # dimension of the observation vector
    x1 <- y[,1]                  # data of channel X
    y1 <- rep(0,ld)
    if( d==2 ) y1 <- y[,2]       # data of channel Y
    if( is.null(lag) )  lag <- as.integer(2*sqrt(ld))    # maximum lag
    lag1 <- lag+1

#  n2p, n : definition    n=2**n2p >= ld+lag1 >= 2**(n2p-1)
    nd <- ld+lag1
    n2p <- 1
    n <- 2**n2p
    while( (n-nd) < 0 ) {
	n2p <- n2p+1
	n <- 2**n2p
    }

    acov <- array(0, dim=c(n,2))      # auto covariance
    ccov21 <- rep(0,n)                # cross covariance
    ccov12 <- rep(0,n)                # cross covariance
    acor <- array(0, dim=c(lag1,2))   # auto correlation
    ccor21 <- rep(0,lag1)             # cross correlation
    ccor12 <- rep(0,lag1)             # cross correlation
    mean <- rep(0,2)                  # mean

    z <- .Call("fftcor",
	as.integer(ld),
	as.integer(lag1),
	as.integer(n),
	as.integer(n2p),
	as.integer(isw),
	as.double(x1),
	as.double(y1))
    
    acv <- array(z[[1L]], dim=c(n,2))
    acov <- array(, dim=c(lag1,2))
    for( i in 1:lag1 ) acov[i,1] <- acv[i,1]
    for( i in 1:lag1 ) acov[i,2] <- acv[i,2]
    
    if( isw == 4 ) {
      fftcor.out  <- list( acov=acov, ccov12=z[[3L]][1:lag1], ccov21=z[[2L]][1:lag1],
	                   acor=array(z[[4L]], dim=c(lag1,2)), ccor12=z[[6L]], ccor21=z[[5L]], mean=z[[7L]] )
      if( plot == TRUE ) {
        par(mfrow=c(2,1))
        plot((0:lag), z[[6L]], type="h", xlab="Lag", ylab="Crosscorrelation ccor12")
        if( lag_axis == TRUE ) abline(h=0, lty=1)
        plot((0:lag), z[[5L]], type="h", xlab="Lag", ylab="Crosscorrelation ccor21")
        if( lag_axis == TRUE ) abline(h=0, lty=1) }
        par(mfrow=c(1,1))
    }

    if( isw != 4 ) fftcor.out  <- list( acov=acov, acor=array(z[[4L]], dim=c(lag1,2)), mean=z[[7L]] )

    return( fftcor.out )
}


mulcor <-
function (y, lag=NULL, plot=TRUE, lag_axis=TRUE)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(lag) )  lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1

    z <- .Call("mulcor",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.integer(lag1))

    cov = array(z[[2L]],dim=c(lag1,d,d))
    cor = array(z[[3L]],dim=c(lag1,d,d))

    mulcor.out <- list( cov=cov, cor=cor, mean=z[[1L]])
    class( mulcor.out ) <- "mulcor"

    if( plot == TRUE ) {
      oldpar <- par(no.readonly=TRUE)
      x <- rep(0,lag1)
      for( i in 1:lag1 ) x[i] <- i-1
      par(mfrow=c(d,d))
      for( j in 1:d )
        for( k in 1:d )  {
          cy <- cor[,j,k]
          plot(x, cy, type="h", xlab="lag", ylab=paste("cor [,",j,",",k,"]"),ylim=c(-1,1))
          if( lag_axis == TRUE ) abline(h=0, lty=1) }
      par(oldpar)
      return( invisible(mulcor.out) )
    } else return( mulcor.out )
}


print.mulcor <- function(x, ...)
{
  lag <- dim(x$cov)[1]-1
  ip <- dim(x$cov)[2]

  cat("\n Mean =" )
  for( i in 1:ip )  cat(sprintf(" %f", x$mean[i]))

  cat("\n\n Autocovariance\n")
  for( i in 1:ip ) {
    cat(sprintf("\n lag\tcov[,%i,%i]\tNormalized\n",i,i))
    for( j in 0:lag ) cat(sprintf(" %i\t%f\t%f\n", j,x$cov[j+1,i,i],x$cor[j+1,i,i]))
  }

  cat("\n\n Cross-covariance\n")
  for( i in 2:ip )
    for( j in 1:(i-1) ) {
      cat(sprintf("\n lag\tcor[,%i,%i]\tNormalized\tcor[,%i,%i]\tNormalized\n",j,i,i,j))
      for( k in 0:lag ) cat(sprintf(" %i\t%f\t%f\t%f\t%f\n", k,x$cov[k+1,j,i],x$cor[k+1,j,i],x$cov[k+1,i,j],x$cor[k+1,i,j]))
    }

}


auspec <-
function (y, lag=NULL, window="Akaike", log=FALSE, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1

    spec1 <- rep(0,lag1)   # spectrum smoothing by window W1
    spec2 <- rep(0,lag1)   # spectrum smoothing by window W2
    stat <- rep(0,lag1)    # test statistics

    z1 <- autcor(y, lag, plot=FALSE)

    z <- .Call("auspec",
	     as.integer(n),
	     as.integer(lag1),
	     as.double(z1$acov))

    if( window == "Akaike" ) spec <- z[[2L]]
    if( window == "Hanning" ) spec <- z[[1L]]

    if( plot == TRUE ) {
      x <- rep(0,lag1)
      for( i in 1:lag1 ) x[i] <- (i-1)/(2*lag)
      if ( log == TRUE ) plot(x, spec, type="l", log="y", xlab="Frequency", ylab="Spectrum" )
      if ( log == FALSE ) plot(x, spec, type="l", xlab="Frequency", ylab="Spectrum")
      auspec.out <- list( stat=z[[3L]] )
    } else auspec.out <- list( spec=spec, stat=z[[3L]] )
    return( auspec.out )
}


mulspe <- 
function (y, lag = NULL, window = "Akaike", plot = TRUE, plot.scale = FALSE) 
{
    n <- nrow(y)
    d <- ncol(y)
    if (is.null(lag)) 
        lag <- as.integer(2 * sqrt(n))
    lag1 <- lag + 1

    z1 <- mulcor(y, lag, plot = FALSE)

    z <- .Call("mulspe", as.integer(n), as.integer(d), as.integer(lag1), 
        as.integer(lag1), as.double(z1$cov))

    if (window == "Akaike") {
        spec <- array(z[[2L]], dim = c(lag1, d, d))
        coh <- array(z[[5L]], dim = c(lag1, d, d))
    }
    if (window == "Hanning") {
        spec <- array(z[[1L]], dim = c(lag1, d, d))
        coh <- array(z[[4L]], dim = c(lag1, d, d))
    }
    cspec <- array(0, dim = c(d, d, lag1))

    mulspe.out <- list(spec = spec, stat = array(z[[3L]], dim = c(lag1, d)), coh = coh)
    class(mulspe.out) <- "mulspe"

    if (plot == TRUE) {
        for (j in 1:d) {
            for (k in 1:j) {
                for (l in 1:lag1) {
                  cspec[j, k, l] <- spec[l, j, k] + (0+1i) * spec[l, k, j]
                  cspec[k, j, l] <- spec[l, j, k] - (0+1i) * spec[l, k, j]
                }
            }
        }
        plot.mulspec(cspec, d, lag, plot.scale)
        return( invisible(mulspe.out) )
    } else   return(mulspe.out)
}


print.mulspe <- function(x, ...)
{
  lag <- dim(x$spec)[1]-1
  ip <- dim(x$spec)[2]

  for( i in 1:ip ) {
    cat(sprintf("\nPower spectrum P(%i,%i)\tSignificance\n",i,i))
    for( j in 0:lag ) cat(sprintf("%i\t%f\t%f\n", j,x$spec[j+1,i,i],x$stat[j+1,i]))
    if( i != 1 ) for( j in 1:(i-1) ) {
      cat(sprintf("\nCross spectrum P(%i,%i)\n",i,j))
      cat("I\tCo-spectrum\tQuad-spectrum\tSimple coherence\n")
      for( k in 0:lag ) cat(sprintf("%i\t%f\t%f\t%f\n", k,x$spec[k+1,i,j],x$spec[k+1,j,i],x$coh[k+1,i,j]))
    }
  }
}


plot.mulspec <- function(spec, d, lag, plot.scale)
{
  par(mfrow = c(d,d))
  lag1 <- lag+1
  x <- rep(0,lag1)
  for (i in 1:lag1) x[i] <- (i-1)/(2*lag)
  dspec <- array(0, dim=c(d,d,lag1))
  for (j in 1:d) for (k in 1:d) {
    if ( j >= k ) {
      dspec[j,k,] = Mod(spec[j,k,])
    } else {
      dspec[j,k,] = Arg(spec[j,k,])
    }
  }

  if(plot.scale == TRUE) {
    mx  <- max(dspec[1,1,])
    for(j in 1:d) for (k in 1:j) if (mx < max(dspec[j,k,])) mx <- max(dspec[j,k,])
    for(j in 1:d) for (k in 1:d)
      if ( j >= k ) {
      ylabs = paste("AmpSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs, ylim=c(0,mx))
    } else {
      ylabs = paste("PhaseSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs, ylim=c(-pi,pi))
    }
  } else {
    for(j in 1:d) for (k in 1:d)
      if ( j >= k ) {
      ylabs = paste("AmpSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs)
    } else {
      ylabs = paste("PhaseSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs)
    }
  }
  par(mfrow=c(1,1))
}


sglfre <-
function (y, lag=NULL, invar, outvar)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))
    lag1 <- lag+1
    z1 <- mulspe(y, lag, window="Akaike", plot=FALSE)

    z <- .Call("sglfre",
	     as.integer(invar),
	     as.integer(outvar),
	     as.integer(n),
	     as.integer(lag1),
	     as.integer(d),
	     as.double(z1$spec))

    sglfre.out <- list( inspec=z[[1L]], outspec=z[[2L]], cspec=z[[3L]], qspec=z[[4L]], gain=z[[5L]],
                          coh=z[[6L]], freqr=z[[7L]], freqi=z[[8L]], errstat=z[[9L]], phase=z[[10L]] )
    return( sglfre.out )
}


mulfrf <-
function (y, lag=NULL, iovar=NULL)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1
    if( is.null(iovar) ) iovar <- c(1:d)
    nv <- length(iovar)-1          # number of input variables
    if( nv == d || nv > d ) stop("number of input variable is smaller than d") 
    for( i in 1: (nv+1) ) {
      if( iovar[i] < 1 ) stop("control variable is greater than or equal to 1")
      if( iovar[i] > d ) stop("control variable is smaller than or equal to d")
    }

    z1 <- mulspe(y, lag, plot=FALSE)

    z <- .Call("mulfrf",
	as.integer(nv),
	as.integer(iovar),
	as.integer(n),
	as.integer(lag1),
	as.integer(d),
	as.double(z1$spec))

    csp=array(z[[1L]], dim=c(d,d,lag1))
    fr=array(z[[2L]], dim=c(nv,lag1))
    fi=array(z[[3L]], dim=c(nv,lag1))
    g=array(z[[4L]], dim=c(nv,lag1))
    ph=array(z[[5L]], dim=c(nv,lag1))
    pch=array(z[[6L]], dim=c(nv,lag1))
    err=array(z[[7L]], dim=c(nv,lag1))
    mch=z[[8L]] 

    mulfrf.out <- list( cospec=csp, freqr=fr, freqi=fi, gain=g, phase=ph, pcoh=pch,
 errstat=err, mcoh=mch )
    return( mulfrf.out )
}


fpeaut <-
function (y, max.order=NULL)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))   # upper limit of model order
    morder <- max.order

    lag <- morder
    lag1 <- lag+1
    z1 <- autcor(y, lag, plot=FALSE)

    z <- .Call("fpeaut",
	as.integer(morder),
	as.integer(n),
	as.double(z1$acov[1]),
	as.double(z1$acov[2:(lag1)]))

    a <- array(z[[11L]],dim=c(morder,morder))
    arcoef <- list()
    for( i in 1:morder ) arcoef[[i]] <- a[1:i,i]

    mo <- z[[9L]]
    fpeaut.out <- list( ordermin=mo, best.ar=arcoef[[mo]], sigma2m=z[[10L]], fpemin=z[[7L]],
                        rfpemin=z[[8L]], ofpe=z[[6L]], arcoef=arcoef, sigma2=z[[1L]], fpe=z[[2L]], rfpe=z[[3L]],
			parcor=z[[4L]], chi2=z[[5L]] )
    return( fpeaut.out )
}


fpec <-
function (y, max.order=NULL, control=NULL, manip=NULL)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)         # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- length(control)       # number of controlled variables
    nman <- length(manip)         # number of manipulated variables
    if( is.null(control) ) {
      ncon <- d
      control <- c(1:d)
      inw <- control
    } else {
      inw <- control
      for( i in 1:ncon ) {
        if( inw[i] < 1 ) stop("control variable is greater than or equal to 1" )
        if( inw[i]  > d ) stop("control variable is smaller than or equal to d" )
      }
    }
    if( is.null(manip) ) {
      nman <- 0
    } else {
      for( i in 1:nman )  {
        if( manip[i] < 1 ) stop("manipulate variable is greater than or equal to 1" )
        if( manip[i]  > d ) stop("manipulate variable is smaller than or equal to d" )
      }
      inw <- c(inw,manip)
    }
    ip <- ncon+nman
    if( ip > d ) stop("length of control and manipulate variables is smaller than or equal to d")


    morder <- max.order
    morder1 <- morder+1
    z1 <- mulcor(y, morder, plot=FALSE)
    cov <- z1$cov[1:morder1,,]

    z <- .Call("fpec7",
	as.integer(n),
	as.integer(morder),
	as.integer(ncon),
	as.integer(ip),
	as.integer(d),
	as.integer(inw),
	as.double(cov))

    cov <- array(z[[1L]], dim=c(morder1,ip,ip))
    cov <- aperm(cov, c(2,3,1))
    arcoef <- array(z[[10L]], dim=c(morder,ncon,ip))
    arcoef <- arcoef[1:z[[5L]],,,drop=F]
    arcoef <- aperm(arcoef, c(2,3,1))

    fpec.out  <- list( cov=cov, fpec=z[[2L]], rfpec=z[[3L]], aic=z[[4L]], ordermin=z[[5L]],
                       fpecmin=z[[6L]], rfpecmin=z[[7L]], aicmin=z[[8L]],
                        perr=array(z[[9L]], dim=c(ncon,ncon)), arcoef=arcoef )
    class( fpec.out ) <- "fpec"
    return( fpec.out )
}

print.fpec <- function(x, ...)
{
  m1 <- dim(x$arcoef)[1]
  m2 <- dim(x$arcoef)[2]
  m3 <- dim(x$arcoef)[3]
  lag <- dim(x$cov)[3]-1

  cat("\n\nCovariance matrix\n")
  print(x$cov)

  cat("m\tFPEC\t\tRFPEC\t\tAIC\n")
  for( i in 0:lag ) cat(sprintf("%i\t%f\t%f\t%f\n",i, x$fpec[i+1],x$rfpec[i+1],x$aic[i+1]))
  cat(sprintf("\n Minimum FPEC = %f\n Minimum RFPEC = %f\n Minimum AIC = %f attained at m = %i\n", x$fpecmin, x$rfpecmin, x$aicmin,  x$ordermin))

  cat("\n\nPrediction error covariance matrix\n")
  print(x$perr)

  cat("\n\nAR coefficient matrix\n\n")
  print(x$arcoef)
}

mulnos <-
function ( y, max.order=NULL, control=NULL, manip=NULL, h)
{
    n <- nrow(y)                    # length of data
    d <- ncol(y)                     # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- length(control)       # number of controlled variables
    nman <- length(manip)         # number of manipulated variables
    n <- nrow(y)        # length of data
    d <- ncol(y)         # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- length(control)       # number of controlled variables
    nman <- length(manip)         # number of manipulated variables
    if( is.null(control) ) {
      ncon <- d
      control <- c(1:d)
      inw <- control
    } else {
      inw <- control
      for( i in 1:ncon ) {
        if( inw[i] < 1 ) stop("control variable is greater than or equal to 1")
        if( inw[i]  > d ) stop("control variable is smaller than or equal to d")
      }
    }
    if( is.null(manip) ) {
      nman <- 0
    } else {
      for( i in 1:nman )  {
        if( manip[i] < 1 ) stop("manipulate variable is greater than or equal to 1")
        if( manip[i]  > d ) stop("manipulate variable is smaller than or equal to d")
      }
      inw <- c(inw,manip)
    }
    ip <- ncon+nman
    if( ip > d ) stop("length of control and manipulate variables is smaller than or equal to d")

    z1 <- fpec(y, max.order, control, manip)
    arcoef <- aperm(z1$arcoef, c(3,1,2))

    z <- .Call("mulnos",
	as.integer(h),
	as.integer(z1$ordermin),
	as.integer(ip),
	as.double(z1$perr),
	as.double(arcoef))

    mulnos.out <- list( nperr=array(z[[1L]], dim=c(ip,ip)), diffr=array(z[[2L]], dim=c(ip,ip,(h+1))),
			integr=array(z[[3L]], dim=c(ip,ip,(h+1))) )
    return( mulnos.out )
}


raspec <-
function (h, var, arcoef=NULL, macoef=NULL, log = FALSE, plot=TRUE)
{
    l <- length(arcoef)
    k <- length(macoef)
    if( is.null(arcoef) )  arcoef <- 0   # coefficient matrix of autoregressive model
    if( is.null(macoef) )  macoef <- 0   # coefficient matrix of moving average model
    macoef <- -macoef

    z <- .Call("raspec",
	as.integer(h),
	as.integer(l),
	as.integer(k),
	as.double(var),
	as.double(arcoef),
	as.double(macoef))

    rspec <- z[[1L]]

    if( plot == TRUE ) {
      x <- rep(0,(h+1))
      for( i in 1:(h+1) ) x[i] <- (i-1)/(2*h)
      if( log == TRUE ) plot(x, rspec, type="l", log="y", xlab="Frequency", ylab="Rational Spectrum")
      if( log == FALSE ) plot(x, rspec, type="l", xlab="Frequency", ylab="Rational Spectrum")
      return( invisible(rspec) )
    } else return( raspec=rspec )
}

mulrsp <-
function (h, d, cov, ar = NULL, ma = NULL, log = FALSE, plot = TRUE, 
    plot.scale = FALSE) 
{
    if (is.null(ar)) {
        l <- 0
        ar <- 0
    }
    else {
        l <- dim(ar)[3]
    }
    if (is.null(ma)) {
        k <- 0
        ma <- 0
    }
    else {
        k <- dim(ma)[3]
    }
    arcoef <- array(0, dim = c(d, d, l))
    macoef <- array(0, dim = c(d, d, k))
    if (l != 0) 
        arcoef <- aperm(ar, c(3, 1, 2))
    if (k != 0) 
        macoef <- aperm(ma, c(3, 1, 2))
    macoef <- -macoef

    z <- .Call( "mulrsp", as.integer(h), as.integer(l), as.integer(d), 
        as.integer(k), as.double(cov), as.double(arcoef), as.double(macoef) )

    mulrsp.out <- list(rspec = array(z[[1L]], dim = c(d,d,(h+1))), scoh = array(z[[2L]], dim = c(d,d,(h+1))))

    if (plot == TRUE) {
        plot.mulspec(mulrsp.out$rspec, d, h, plot.scale)
        return( invisible(mulrsp.out) )
    } else return(mulrsp.out)
}


optdes <-
function (y, max.order=NULL, ns, q, r)
{
    n <- nrow(y)       # length of data
    d <- ncol(y)       # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- dim(q)[1]  # number of controlled variables
    nman <- dim(r)[1]  # number of manipulated variables
    control <- c(1:ncon)
    manip <- c((ncon+1):(ncon+nman))
    z1 <- fpec(y, max.order, control, manip)
    order <- z1$ordermin
    ao <- z1$arcoef
    ao <- aperm(ao, c(3,1,2))
    ao <- aperm(ao, c(2,1,3))
    dim(ao) <- c(ncon*order,(ncon+nman))
    a <- ao[,(1:ncon)]
    b <- ao[,(ncon+1):(ncon+nman)]

    z <- .Call("optdes",
	as.integer(ncon),
	as.integer(nman),
	as.integer(ns),
	as.integer(order),
	as.double(q),
	as.double(r),
	as.double(z1$perr),
	as.double(a),
	as.double(b))

    optdes.out <- list( perr=z1$perr, trans=a, gamma=b, gain=array(z[[1L]], dim=c(nman,order*ncon)) )
    return( optdes.out )
}


optsim <-
function (y, max.order=NULL, ns, q, r, noise=NULL, len, plot=TRUE)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ir <- dim(q)[1]     # number of controlled variables
    il <- dim(r)[1]     # number of manipulated variables

    z1 <- optdes(y, max.order, ns, q, r)
    trans <- z1$trans
    gamma <- z1$gamma
    gain <- z1$gain
    order <- dim(gain)[2]/ir
    if( is.null(noise) ) noise <- wnoise(len, z1$perr, plot=FALSE)     # white noise

    z <- .Call("optsim",
	as.integer(ns),
	as.integer(order),
	as.integer(ir),
	as.integer(il),
	as.double(trans),
	as.double(gamma),
	as.double(gain),
	as.double(noise))

    convar <- array(z[[1L]],dim=c(ir,ns))
    manvar <- array(z[[2L]],dim=c(il,ns))

    if( plot == TRUE ) {
      nc <- max(ir,il)
      par(mfcol=c(2,nc))
      for( i in 1:ir ) {
        plot(convar[i,], type="h", xlab="Step",ylab=paste("Controlled variables X(",i,")"))
        abline(h=0, lty=1) }
      for( i in 1:il ) {
        plot(manvar[i,], type="h", xlab="Step",ylab=paste("Manipulated variables Y(",i,")"))
        abline(h=0, lty=1) }
      par(mfrow=c(1,1))
      optsim.out <- list( trans=trans, gamma=gamma, gain=gain, 
			xmean=z[[3L]], ymean=z[[4L]], xvar=z[[9L]], yvar=z[[10L]],
			x2sum=z[[5L]], y2sum=z[[6L]], x2mean=z[[7L]], y2mean=z[[8L]] )
    } else {
      optsim.out <- list( trans=trans, gamma=gamma, gain=gain, convar=convar, manvar=manvar,
			xmean=z[[3L]], ymean=z[[4L]], xvar=z[[9L]], yvar=z[[10L]],
			x2sum=z[[5L]], y2sum=z[[6L]], x2mean=z[[7L]], y2mean=z[[8L]] )
    }
    return( optsim.out )
}


wnoise <-
function (len, perr, plot=TRUE)
{
    ir <- 1                           # number of controlled variables
    if( is.array(perr) ) ir <- ncol(perr)

    z <- .Call("wnoise",
	as.integer(len),
	as.integer(ir),
	as.double(perr))

    if( ir != 1 ) wnoise <- array(z[[1L]],dim=c(ir,len))
    if( ir == 1 ) wnoise <- z[[1L]][1:len]

    if( plot == TRUE ) {
#      nc <- as.integer((ir+1)/2)
#      if( ir != 1 ) par(mfrow=c(2,nc))
      par(mfcol=c(ir,1))
      if( ir == 1 ) {
        plot(wnoise, type="h", main="White Noise", xlab="n", ylab=expression(y[n]))
        abline(h=0, lty=1)
      } else {
        plot(wnoise[1,], type="h", main=expression(paste("White noise ", y[n])), xlab="n", ylab="d=1")
        abline(h=0, lty=1)
        for( i in 2:ir ) {
          plot(wnoise[i,], type="h", xlab="n", ylab=paste("d=",i))
          abline(h=0, lty=1)
        }
      }
      par(mfrow=c(1,1))
    } else  return( wnoise=wnoise )
}

#####   TIMSAC74   #####

autoarmafit <-
function (y, max.order=NULL)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))

    morder <- max.order
    autcv <- autcor(y, morder, plot=FALSE)$acov   # covariance sequence
    z1 <- canarm(y, morder, plot=FALSE)
    inc <- 1               # total number of case
    arcoef1 <- z1$arcoef   # initial estimates of AR-coefficients
    macoef1 <- z1$macoef   # initial estimates of MA-coefficients
    p1 <- length(arcoef1)  # initial AR order
    q1 <- length(macoef1)  # initial MA order

    icst <- 190
    mmax <- morder
    lmax <- morder+icst+mmax
    nmax <- 25

    z <- .Call("autarm",
	as.integer(n),
	as.integer(morder+1),
	as.double(autcv),
	as.integer(inc),
	as.integer(p1),
	as.double(arcoef1),
	as.integer(q1),
	as.double(macoef1),
	as.integer(lmax),
	as.integer(mmax),
	as.integer(nmax) )

    nc <- z[[1L]]
    p <- z[[2L]][1:nc]
    q <- z[[4L]][1:nc]
    a <- array(z[[3L]], dim=c(mmax,nmax))
    b <- array(z[[5L]], dim=c(mmax,nmax))
    std <- array(z[[6L]], dim=c(mmax,nmax))
    v <- rep(0,nc)
    gr <- array(z[[8L]], dim=c(mmax,nmax))
    aaic <- z[[9L]][1:nc]
    arma <- vector("list", nc)
    arstd <- list()
    mastd <- list()
    grad <- list()

    aicorder <- sort(aaic, index.return=TRUE)

    for ( i in 1:nc ) {
        j <- aicorder$ix[i]
        if( p[j] == 0 ) {
          ar <- NULL
        } else {
          ar <- -a[(1:p[j]),j]
        }
        if( q[j] == 0 ) {
          ma <- NULL
        } else {
          ma <- -b[(1:q[j]),j]
        }
        arma[[i]] <- list(arcoef=ar, macoef=ma)
        arstd[[i]] <- std[(q[j]+1):(p[j]+q[j]),j]
        mastd[[i]] <- std[1:q[j],j]
        v[i] <- z[[7L]][j]
        grad[[i]] <- gr[(1:(p[j]+q[j])),j]
    }

    best.model <- list( arcoef=arma[[1]]$arcoef, macoef=arma[[1]]$macoef, arorder=z[[11L]], maorder=z[[12L]] )

    model <- list()
    for ( i in 1:nc )
      model[[i]] <- list( arcoef=arma[[i]]$arcoef, macoef=arma[[i]]$macoef, arstd=arstd[[i]], mastd=mastd[[i]], v=v[i],
                 aic=aicorder$x[i], grad=grad[[i]] )

    autoarmafit.out <- list( best.model=best.model, model=model )
    class( autoarmafit.out ) <- "autoarmafit"
    return( autoarmafit.out )
}

print.autoarmafit <- function(x, ...)
{
  nc <- length(x$model)
  z <- x$model

  for( i in 1:nc) {
    iq <- length(z[[i]]$arcoef)
    ip <- length(z[[i]]$macoef)

    cat(sprintf("\n\nCase No. %i\n", i))
    if( iq != 0 ) {
      cat("\n AR coefficient\tStandard deviation\n")
      for( j in 1:iq ) cat(sprintf(" %f\t%f\n", z[[i]]$arcoef[j], z[[i]]$arstd[j]))
    } 
    if( ip != 0 ) {  
      cat("\n MA coefficient\tStandard deviation\n")
      for( j in 1:ip ) cat(sprintf(" %f\t%f\n", z[[i]]$macoef[j], z[[i]]$mastd[j]))
    }
    cat(sprintf("\n AIC\t%f\n", z[[i]]$aic))
    cat(sprintf(" Innovation variance\t%f\n",z[[i]]$v))
    cat(" Final gradient")
    for( j in 1:(iq+ip) ) cat(sprintf("\t%e", z[[i]]$grad[j]))
    cat("\n")
  }

  zz <- x$best.model
  cat("\n\nBest ARMA model")
  cat(sprintf("\n AR coefficient (order = %i)", zz$arorder))
  for( j in 1:zz$arorder ) cat(sprintf("\t%f", zz$arcoef[j]))
  cat(sprintf("\n MA coefficient (order = %i)", zz$maorder))
  for( j in 1:zz$maorder ) cat(sprintf("\t%f", zz$macoef[j]))
  cat("\n\n")
}

armafit <-
function (y, model.order)
{
    n <- length(y)
    lag <- as.integer(2*sqrt(n))     # maximum lag
    lag1 <- lag+1

    inc <- 1             # total number of case
    p1 <- rep(0,2)
    q1 <- rep(0,2)

    autcv <- autcor(y, lag, plot=FALSE)$acov		# covariance sequence
    z1 <- canarm(y, lag, plot=FALSE)
    arcoef1 <- z1$arcoef	# initial estimates of AR-coefficients
    macoef1 <- z1$macoef	# initial estimates of MA-coefficients
    p1[1] <- length(arcoef1)	# initial AR order
    q1[1] <- length(macoef1)	# initial MA order

    if( length(model.order) == 2 ) {
      inc <- inc+1
      p1[2] <- model.order[1]     # AR order to be fitted successively
      q1[2] <- model.order[2]     # MA order to be fitted successively
    }

    icst <- 190
    mmax <- 50
    lmax <- lag+icst+mmax
    nmax <- 25	         # the limit of the total cases

    z <- .Call("autarm",
	as.integer(n),
	as.integer(lag1),
	as.double(autcv),
	as.integer(inc),
	as.integer(p1),
	as.double(arcoef1),
	as.integer(q1),
	as.double(macoef1),
	as.integer(lmax),
	as.integer(mmax),
	as.integer(nmax) )

    nc <- z[[1L]]
    p <- z[[2L]][1:nc]
    q <- z[[4L]][1:nc]
    a <- array(z[[3L]], dim=c(mmax,nmax))
    b <- array(z[[5L]], dim=c(mmax,nmax))
    std <- array(z[[6L]], dim=c(mmax,nmax))
    gr <- array(z[[8L]], dim=c(mmax,nmax))
    aaic <- z[[9L]][1:nc]
    aic <- NULL
    arcoef <- NULL
    macoef <- NULL
    arstd <- NULL
    mastd <- NULL
    grad <- NULL

    aicorder <- sort(aaic, index.return=TRUE)
    order.maice <- aicorder$ix[1]

    arorder <- model.order[1]
    maorder <- model.order[2]

    for( i in 1:nc ) {
        if ( is.null(aic) == FALSE ) break
        j <- aicorder$ix[i]
        if( p[j]==arorder  &&  q[j]==maorder) {
          arcoef <- -a[(1:p[j]),j]
          macoef <- -b[(1:q[j]),j]
          arstd <- std[(q[j]+1):(p[j]+q[j]),j]
          mastd <- std[1:q[j],j]
          v <- z[[7L]][j]
          aic <- z[[9L]][j]
          grad <- gr[(1:(p[j]+q[j])),j]
        }
    }

    armafit.out <- list( arcoef=arcoef, macoef=macoef, arstd=arstd, mastd=mastd, v=v, aic=aic, grad=grad )
    return( armafit.out )
}


bispec <-
function (y, lag=NULL, window="Akaike", log=FALSE, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1

    z1 <- thirmo(y, lag, plot=FALSE)
    cv <- z1$acov			# autocovariance

    tmnt <- array(0,dim=c(lag1,lag1))
    for( i in 1:lag1 ) tmnt[i,1:i] <- z1$tmomnt[[i]]   # third order moments

    z <- .Call("bispec",
	as.integer(n),
	as.integer(lag),
	as.double(cv),
	as.double(tmnt))

    if( window == "Akaike" ) pspec <- z[[2L]]
    if( window == "Hanning" ) pspec <- z[[1L]]

    if( plot == TRUE ) {
      x <- rep(0,lag1)
      for( i in 1:lag1 ) x[i] <- (i-1)/(2*lag)
      ylab = paste("Spectrum smoothing by ",window," window")
      if( log == TRUE ) plot(x, pspec, type="l", log="y", ylab=ylab, xlab="Frequency")
      if( log == FALSE ) plot(x, pspec, type="l", ylab=ylab, xlab="Frequency") }

    bispec.out <- list( pspec=pspec, sig=z[[3L]], cohe=array(z[[4L]], dim=c(lag1,lag1)),
			breal=array(z[[5L]], dim=c(lag1,lag1)), bimag=array(z[[6L]], dim=c(lag1,lag1)), exval=z[[7L]] )
    return( bispec.out )
}


canarm <-
function (y, lag=NULL, max.order=NULL, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))
    morder <- lag

    z2 <- autcor(y, morder, plot=FALSE)
    autcv <- z2$acov

    if( is.null(max.order) ) max.order <- morder
    mmax <- max.order
#    ifpl <- as.integer(3.0*sqrt(n))
    ifpl <- min(mmax, morder)
    l1 <- ifpl+1			# upper limit of the model order +1
    nmax <- 2*morder+1

    z1 <- .Call("canarm",
	as.integer(n),
	as.integer(morder+1),
	as.double(autcv),
	as.integer(l1),
	as.integer(mmax),
	as.integer(nmax))

    ar <- z1[[1L]]
    mo <- z1[[5L]]
    parcor <- z1[[6L]][1:(l1-1)]
    nc <- z1[[7L]]
    m1 <- z1[[8L]][1:nc]
    m2 <- z1[[9L]][1:nc]
    w <- array(z1[[10L]], dim=c(mmax,mmax,nc))
    z <- array(z1[[11L]], dim=c(mmax,nc))
    Rs <- array(z1[[12L]], dim=c(mmax,nc))
    chi <- array(z1[[13L]], dim=c(mmax,nc))
    ndt <- array(z1[[14L]], dim=c(mmax,nc))
    dicp <- array(z1[[15L]], dim=c(mmax,nc))
    cweight <- list()
    cR <- list()
    Rsquar <- list()
    chisquar <- list()
    ndf <- list()
    dic <- list()
    for( i in 1:nc )  cweight[[i]] <- w[(1:m1[i]),(1:m1[i]),i]
    for( i in 1:nc )  cR[[i]] <- z[(1:m1[i]),i]
    for( i in 1:nc )  Rsquar[[i]] <- Rs[(1:m1[i]),i]
    for( i in 1:nc )  chisquar[[i]] <- chi[(1:m1[i]),i]
    for( i in 1:nc )  ndf[[i]] <- ndt[(1:m1[i]),i]
    for( i in 1:nc )  dic[[i]] <- dicp[(1:m1[i]),i]
    k <- z1[[18L]]
    l <- z1[[20L]]

    if( plot == TRUE ) {
      plot(parcor, type="h", ylab="Partial Autocorrelation", xlab="Lag")
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n), lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3) }

    canarm.out <- list( arinit=-ar[1:mo], v=z1[[2L]][1:l1], aic=z1[[3L]][1:l1], aicmin=z1[[4L]],
                        order.maice=mo, parcor=parcor, nc=nc, future=m1, past=m2, cweight=cweight,
                        canocoef=cR, canocoef2=Rsquar, chisquar=chisquar, ndf=ndf, dic=dic, dicmin=z1[[16L]][1:nc],
                        order.dicmin=z1[[17L]][1:nc], arcoef=-z1[[19L]][1:k], macoef=-z1[[21L]][1:l] )
    return( canarm.out )
}


canoca <-
function (y)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of Y(I)
    lag <- as.integer(2*sqrt(n))     # maximum lag
    lag1 <- lag+1
    inw <- c(1:d)       # inw(k)=j means that the k-th component of y(i) is the j-th component of the original record z(i)

    z2 <- mulcor(y, lag, plot=FALSE)
    cov <- z2$cov       # covariance matrix

    lmax <- 12
    mj0 <- lmax+1
    mj1 <- mj0*d

    l <- 0		# upper bound of AR-order
    aic <- rep(0,mj0)	# AIC
    oaic <- 0		# minimum AIC
    mo <- 0		# MAICE AR-model order                 
    v <- array(0, dim=c(d,d))	# innovation variance
    ac <- array(0, dim=c(mj0,d,d))	# autoregressive coefficients
    nc <- 0				# number of cases
    m1 <- rep(0,mj1)			# number of variable in the future set
    m2 <- rep(0,mj1)			# number of variables in the past set
    w <- array(0,dim=c(mj1,mj1,mj1))	# future set canonical weight
    z <- array(0,dim=c(mj1,mj1))	# canonical R
    Rs <- array(0,dim=c(mj1,mj1))	# R-squared
    chi <- array(0,dim=c(mj1,mj1))	# chi-square
    ndt <- array(0,dim=c(mj1,mj1))	# N.D.F
    dic <- array(0,dim=c(mj1,mj1))	# DIC(=CHI**2-2*D.F)
    dicm <- rep(0,mj1)		# minimum DIC
    po <- rep(0,mj1)			# order of minimum DIC
    f <- array(0,dim=c(mj1,mj1))	# transition matrix F
    k <- 0				# number of structual characteristic vector
    nh <- rep(0,mj1)			# structual characteristic vector
    g <- array(0, dim=c(mj1,d))	# input matrix
    ivf <- 0				# number of vector vf
    vf <- rep(0,mj1*mj1)		# F matrix in vector form

    z1 <- .Call("canoca",
	as.integer(d),
	as.integer(inw),
	as.integer(n),
	as.integer(lag1),
	as.integer(d),
	as.double(cov),
	as.integer(lmax),
	as.integer(mj0),
	as.integer(mj1) )

    l <- z1[[1L]]
    mo <- z1[[4L]]
    ac <- array(z1[[6L]], dim=c(mj0,d,d))
    arcoef <- array(,dim=c(d,d,mo))
    for( i in 1:mo ) arcoef[,,i] <- -ac[i,(1:d),(1:d)]
    nc <- z1[[7L]]
    m1 <- z1[[8L]][1:nc]
    m2 <- z1[[9L]][1:nc]
    w <- array(z1[[10L]], dim=c(mj1,mj1,mj1))
    z <- array(z1[[11L]], dim=c(mj1,mj1))
    Rs <- array(z1[[12L]], dim=c(mj1,mj1))
    chi <- array(z1[[13L]], dim=c(mj1,mj1))
    ndt <- array(z1[[14L]], dim=c(mj1,mj1))
    dicp <- array(z1[[15L]], dim=c(mj1,mj1))
    f <- array(z1[[18L]], dim=c(mj1,mj1))
    cweight <- list()
    cR <- list()
    Rsquar <- list()
    chisquar <- list()
    ndf <- list()
    dic <- list()
    matF <- list()
    for( i in 1:nc )  cweight[[i]] <- w[(1:m1[i]),(1:m1[i]),i]
    for( i in 1:nc )  cR[[i]] <- z[(1:m1[i]),i]
    for( i in 1:nc )  Rsquar[[i]] <- Rs[(1:m1[i]),i]
    for( i in 1:nc )  chisquar[[i]] <- chi[(1:m1[i]),i]
    for( i in 1:nc )  ndf[[i]] <- ndt[(1:m1[i]),i]
    for( i in 1:nc )  dic[[i]] <- dicp[(1:m1[i]),i]
    for( i in 1:nc )  matF[[i]] <- f[1:(m1[i]-1),i]
    k <- z1[[19L]]
    g <- array(z1[[21L]], dim=c(mj1,d))
    ivf <- z1[[22L]]

    canoca.out <- list( aic=z1[[2L]][1:(l+1)], aicmin=z1[[3L]], order.maice=mo, v=array(z1[[5L]], dim=c(d,d)),
			arcoef=arcoef, nc=nc, future=m1, past=m2, cweight=cweight, canocoef=cR, canocoef2=Rsquar,
			chisquar=chisquar, ndf=ndf, dic=dic, dicmin=z1[[16L]][1:nc], order.dicmin=z1[[17L]][1:nc],
                        matF=matF, vectH=z1[[20L]][1:k], matG=g[(1:k),(1:d)], vectF=z1[[23L]][1:ivf] )
    return( canoca.out )
}


covgen <-
function (lag, f, gain, plot=TRUE)
{
    k <- length(f)	# number of data points

    z <- .Call("covgen",
	as.integer(lag),
	as.integer(k),
	as.double(f),
	as.double(gain))

    if( plot == TRUE ) {
      plot((0:lag), z[[2L]], type="h", ylab="Auto Covariance Normalized", xlab="Lag")
      abline(h=0, lty=1) }

    covgen.out <- list( acov=z[[1L]], acor=z[[2L]] )
    return( covgen.out )
}


markov <-
function (y)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1
    icont <- 2		# output control (0:ARMA coefficients, 1:for SIMCON input, 2: for both)

    z1 <- mulcor(y, lag, plot=FALSE)
    cov <- z1$cov		# covariance matrix

    z2 <- canoca(y)
    nh <- z2$vectH		# structural characteristic vector
    k <- length(nh)		# dimension of the state vector
    vectF <- z2$vectF		# initial estimate of the vector of free parameters in F
    nvf <- length(vectF)		# length of vector vf
    matGi <- rep(z2$matG,1)     # initial estimates of the free parameters in G

    mj3 <- max(lag1,100)
    mj4 <- nvf+k*d
    mj6 <- 2*(k+d)-1
    mj7 <- as.integer((k-1)/d+1)

    z <- .Call("markov",
	as.integer(n),
	as.integer(lag1),
	as.integer(d),
	as.double(cov),
	as.integer(k),
	as.integer(nh),
	as.integer(nvf),
	as.double(vectF),
	as.double(matGi),
	as.integer(icont),
	as.integer(mj3),
	as.integer(mj4),
	as.integer(mj6),
	as.integer(mj7) )

    ngr <- z[[5L]]
    vd <- array(z[[10L]], dim=c(mj4,mj4))
    iqm <- z[[11L]]
    bm <- array(z[[12L]], dim=c(d,d,mj7))
    au <- array(z[[13L]], dim=c(d,d,mj7))
    zz <- array(z[[14L]], dim=c(d,d,mj7))

    arcoef <- array(-bm, dim=c(d,d,iqm))
    macoef <- array(-zz, dim=c(d,d,iqm-1))

    markov.out <- list( id=z[[1L]], ir=z[[2L]], ij=z[[3L]], ik=z[[4L]], grad=z[[6L]][1:ngr], matFi=array(z[[7L]], dim=c(k,k)),
			matF=array(z[[8L]], dim=c(k,k)), matG=array(z[[9L]], dim=c(k,d)), davvar=vd[(1:ngr),(1:ngr)],
			arcoef=arcoef, impulse=array(au, dim=c(d,d,iqm-1)),
			macoef=macoef, v=array(z[[15L]], dim=c(d,d)), aic=z[[16L]] )
    return( markov.out )
}


nonst <-
function (y, span, max.order=NULL, plot=TRUE)
{
    n <- length(y)	# length of data
    if( span < 1 ) span <- n
    ns <- as.integer(n/span)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # highest order of AR model
    morder <- max.order

    z <- .Call("nonst",
	as.integer(n),
	as.integer(span),
	as.double(y),
	as.integer(ns),
	as.integer(morder))

    coef <- array(z[[2L]], dim=c(morder,ns))
    arcoef <- list() 
    for (i in 1:ns )  arcoef[[i]] <- coef[1:z[[1L]][i],i]
    pspec <- array(z[[9L]], dim=c(121,ns))

    if( plot == TRUE ) {
      oldpar <- par(no.readonly=TRUE)
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      par(mfrow=c(ns,1))
      for( i in 1:ns ) {
        plot(x, pspec[,i], type="l", main=paste("y(", z[[7L]][i], "),...,y(", z[[8L]][i], ")"),
             xlab="Frequency", ylab="Power Spectrum") }
      par(oldpar)
      nonst.out <- list( ns=ns, arcoef=arcoef, v=z[[3L]], aic=z[[4L]], daic21=z[[5L]], daic=z[[6L]],
		       init=z[[7L]], end=z[[8L]] )
    } else {
      nonst.out <- list( ns=ns, arcoef=arcoef, v=z[[3L]], aic=z[[4L]], daic21=z[[5L]], daic=z[[6L]],
		       init=z[[7L]], end=z[[8L]], pspec=pspec )
    }
    class( nonst.out ) <- "nonst"
    return( nonst.out )
}


print.nonst <- function(x, ...)
{
  for( i in 1:x$ns ) {
    cat("\n\n..........  Current model  ..........\n\n")
    cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n", x$init[i],x$end[i]))
    cat(sprintf(" Innovation variance = %e\n", x$v[i]))
    cat(sprintf(" aic = %f\n", x$aic[i]))
    if( i > 1 ) {
      cat(sprintf(" aic2-aic1 = %f\n", x$daic21[i]))
      nd <- x$end[i]-x$init[i]+1
      cat(sprintf(" daic/%i = %f\n", nd,x$daic[i]))
    }
    mf <- length(x$arcoef[[i]])
    cat(sprintf("\n Autoregressive coefficient ( order %i )\n", mf))
    print(x$arcoef[[i]])
  }
}


prdctr <-
function (y, r, s, h, arcoef, macoef=NULL, impulse=NULL, v, plot=TRUE)
{
    if (is.array(y)) {
      n <- nrow(y)      # length of data
      d <- ncol(y)	# dimension of vector y(i)
    } else {
      n <- length(y)
      d <- 1
      y <- array(y, dim=c(n,1))
    }

    if( is.null(arcoef) ) stop( "'arcoef' must be numeric" )
    if(is.array(arcoef)) {	# AR-coefficient matrices
        p <- length(arcoef)[3]
        arcoef <- -arcoef
    } else {
        p <- length(arcoef)
        arcoef <- array(-arcoef, dim=c(1,1,p))
    }

    if( is.null(macoef) && is.null(impulse) ) stop( "'macoef' or ' impulse' must be numeric" )
    jsw <- 0
    if( is.null(macoef) ) {
      jsw <- 1
      if( is.array(impulse) ) {	# impulse response matrices
        q <- dim(impulse)[3]
      } else {
        q <- length(impulse)
        impulse <- array(impulse,dim=c(1,1,q))
      }
      macoef <- array(0, dim=c(d,d,q))
    } else {
      if( is.array(macoef) ) {	# MA-coefficient matrices
        q <- dim(macoef)[3]
        macoef <- -macoef
      } else {
        q <- length(macoef)
        macoef <- array(-macoef,dim=c(1,1,q))
      }
      if( is.null(impulse) ) impulse <- array(0, dim=c(d,d,q))
    }

    z <- .Call("prdctr",
	as.integer(n),
	as.integer(r),
	as.integer(s),
	as.integer(h),
	as.integer(d),
	as.integer(p),
	as.integer(q),
	as.integer(jsw),
	as.double(y),
	as.double(arcoef),
	as.double(macoef),
	as.double(impulse),
	as.double(v))

#    yreal <- array(z$yreal, dim=c(s+h,d))
#    yori <- array(z$yori, dim=c(h+1,d))
#    for(i in 1:(n-s+1)) yreal[s+i-1,] <- yori[i,]
#    for(j in 1:d) yreal[,j] <- yreal[(1:n),j]

    predct <- array(z[[3L]], dim=c(s+h,d))
    for(i in 1:(r-1)) predct[i,] <- NA

#    ys <- array(z$ys, dim=c(n,d))
#    for(i in 1:(r-1)) ys[i,] <- NA
#    for(i in s:n) ys[i,] <- NA
    ys <- array(NA, dim=c(n,d))
    for(j in 1:d)
      for(i in r:n) ys[i,j] <- y[i,j]-predct[i,j]

    pstd <- array(z[[5L]], dim=c(s+h,d))
    pstd[1:s-1,] <- NA
    p2std <- array(z[[6L]], dim=c(s+h,d))
    p2std[1:s-1,] <- NA
    p3std <- array(z[[7L]], dim=c(s+h,d))
    p3std[1:s-1,] <- NA
    mstd <- array(z[[8L]], dim=c(s+h,d))
    mstd[1:s-1,] <- NA
    m2std <- array(z[[9L]], dim=c(s+h,d))
    m2std[1:s-1,] <- NA
    m3std <- array(z[[10L]], dim=c(s+h,d))
    m3std[1:s-1,] <- NA

    if( plot == TRUE ) {
      par(mfrow=c(1,d))
      for( i in 1:d ) {
        ymin <- min(y[(1:n),i], predct[r:(s+h),i])
        ymax <- max(y[(1:n),i], predct[r:(s+h),i])
        plot(y[,i], type="l", xlim=c(0,(s+h)), ylim=c(ymin,ymax), xlab="n", main="real data & predicted values")
        par(new=TRUE)
        plot(predct[,i], type="l", xlim=c(0,(s+h)), ylim=c(ymin,ymax), col="red", xlab="", ylab="") }
      par(mfrow=c(1,1))
    }

    prdctr.out <- list( predct=predct, ys=ys, pstd=pstd, p2std=p2std, p3std=p3std,
                        mstd=mstd, m2std=m2std, m3std=m3std )

    class( prdctr.out ) <- "prdctr"
    return( prdctr.out )
}

print.prdctr <- function(x, ...)
{
  n <- dim(x$ys)[1]
  d <- dim(x$ys)[2]
  sh <- dim(x$predct)[1]
  for (i in n:1 ) {
    if ( is.na(x$ys[i,1]) == FALSE ) r <- i
    if ( is.na(x$pstd[i,1]) == FALSE ) s <- i
  }

  for ( id in 1:d ) {
    cat(sprintf("\n\n d = %i\n\n", id))
    cat(" n\tpredct\t\tys\t\tpstd\t\tp2std\t\tp3std\n")
    cat("\t\t\t\t\tmstd\t\tm2std\t\tm3std\n")
    for( i in r:(s-1) ) cat(sprintf(" %i\t%f\t%f\n", i,x$predct[i,id],x$ys[i,id]))
    for( i in s:n ) {
      cat(sprintf(" %i\t%f\t%f\t%f\t%f\t%f\n", i,x$predct[i,id],x$ys[i,id],x$pstd[i,id],x$p2std[i,id],x$p3std[i,id]))
      cat(sprintf("\t\t\t\t\t%f\t%f\t%f\n", x$mstd[i,id],x$m2std[i,id],x$m3std[i,id])) }
    for( i in (n+1):sh ) {
      cat(sprintf(" %i\t\t\t%f\t%f\t%f\t%f\n", i,x$pstd[i,id],x$predct[i,id],x$p2std[i,id],x$p3std[i,id]))
      cat(sprintf("\t\t\t\t\t%f\t%f\t%f\n", x$mstd[i,id],x$m2std[i,id],x$m3std[i,id])) }
  }

}


simcon <-
function (span, len, r, arcoef, impulse, v, weight)
{
    arcoef <- -arcoef	# matrices of autoregressive coefficients
    d <- dim(arcoef)[1]		# dimension of Y(I)
    k <- dim(arcoef)[3]		# order of the process

    z <- .Call("simcon",
	as.integer(d),
	as.integer(k),
	as.integer(span),
	as.integer(len),
	as.integer(r),
	as.double(arcoef),
	as.double(impulse),
	as.double(v),
	as.double(weight))

    simcon.out <- list( gain=array(z[[3L]], dim=c(r,k*d)), ave=z[[4L]][1:d], var=z[[5L]][1:d], std=z[[6L]][1:d],
			bc=array(z[[1L]], dim=c(k*d,r)), bd=array(z[[2L]], dim=c(k*d,d-r)) )
    return( simcon.out )
}


thirmo <-
function (y, lag=NULL, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1

    z <- .Call("thirmo",
	as.integer(n),
	as.integer(lag),
	as.double(y))

    mnt <- array(z[[4L]], dim=c(lag1,lag1))
    tmomnt <- list()
    for( i in 1:lag1 ) tmomnt[[i]] <- mnt[i,1:i]

    if( plot == TRUE ) {
      plot((0:lag), z[[3L]], type="h", ylab="Auto Covariance Normalized", xlab="Lag")
      abline(h=0, lty=1)
    }

    thirmo.out <- list( mean=z[[1L]], acov=z[[2L]], acor=z[[3L]], tmomnt=tmomnt )
    return( thirmo.out )
}


#####   TIMSAC78   #####

blocar <-
function (y, max.order=NULL, span, plot=TRUE)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of the order of AR model

    morder <- max.order
    if( span < 0 ) span <- n
    ns <- as.integer((n-morder+span-1)/span)

    z <- .Call("blocar",
	as.double(y),
	as.integer(n),
	as.integer(morder),
	as.integer(span),
	as.integer(ns))

    aic <- list()
    bweight <- list()
    pacoef <- list()
    arcoef <- list()
    aic[[1]] <- NA
    bweight[[1]] <- NA
    for( i in 2:ns ) {
        j <- ((i-1)*ns+1):((i-1)*ns+i)
        aic[[i]] <- z[[3L]][j]
        bweight[[i]] <- z[[4L]][j]
    }
    for( i in 1:ns ) pacoef[[i]] <- z[[5L]][((i-1)*morder+1):(i*morder)]
    for( i in 1:ns ) arcoef[[i]] <- z[[6L]][((i-1)*morder+1):(i*morder)]
    pspec <- array(z[[10L]], dim=c(121,ns))

    if( plot == TRUE ) {
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      ymin <- min(pspec)
      ymax <- max(pspec)
      par(mfrow=c(ns,1))
      nk <- ns
      if( ns == 4 ) {
        par(mfrow=c(2,2))
        nk <- 4
      }
      if( ns > 4 ) {
        par(mfrow=c(3,2))
        nk <- 6
      }
      if( ns > 6 ) {
        par(mfrow=c(3,3))
        nk <- 9
      }
      nn <- 0
      for(i in 1:ns) {
        if( nn == nk ) {
          par(ask=TRUE)
          nn <- 0
        }
        plot(x, pspec[,i], type="l", ylim=c(ymin,ymax), main=paste("y(", z[[8L]][i], "),...,y(", z[[9L]][i], ")"),
                          xlab="Frequency", ylab="Power Spectrum")
        nn <- nn+1
      }
      par(mfrow=c(1,1))
    }


    blocar.out <- list( mean=z[[1L]], var=z[[2L]], aic=aic, bweight=bweight, pacoef=pacoef, arcoef=arcoef,
			v=z[[7L]], init=z[[8L]], end=z[[9L]], pspec=pspec )
    class( blocar.out ) <- "blocar"
    return( blocar.out )
}


print.blocar <- function(x, ...)
{
  cat(sprintf("\n Mean\t%f\n", x$mean))
  cat(sprintf(" Variance\t%f\n", x$var))

  ns <- length(x$bweight)
  for( i in 1:ns ) {
    if( i != 1 ) {
      cat("\n\nAR-Model fitted to\t! Bayesian weights\t! AIC with respect to the present data\n")
      cat("--------------------------------------------------------------------------------------\n")
      cat(sprintf("Current block\t\t! %f\t\t! %f\n", x$bweight[[i]][1], x$aic[[i]][1]))
      for( k in 1:(i-1) ) cat(sprintf("%i period former block\t! %f\t\t! %f\n", k,x$bweight[[i]][k+1],x$aic[[i]][k+1]))
    }

    morder <- length(x$arcoef[[i]])
    cat("\n\n..........  Current model  (Average by the Bayesian weight)   ..........\n\n")
    cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n", x$init[i],x$end[i]))
    cat(sprintf(" Innovation variance = %e\n\n", x$v[i]))
    cat(sprintf(" AR coefficients ( order %i ) \n", morder ))
    print(x$arcoef[[i]])
  }
  cat("\n")
}


blomar <-
function (y, max.order=NULL, span)
{max.order
    n <- nrow(y)
    d <- ncol(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of the order of AR-model
    icflag <- 0
    if( max.order > n/(2*d) )  icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))
    morder <- max.order

    calb<-rep(1,d)   # calibration constant for channel j (j=1,d)
    if( span < 1 ) span <- n
    ns <- as.integer((n-morder+span-1)/span)

    z <- .Call("blomar",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.double(calb),
	as.integer(morder),
	as.integer(span),
	as.integer(ns))

    nns <- z[[10L]]
    bw <- array(z[[3L]], dim=c(ns,ns))
    bweight <- list()
    bweight[[1]] <- NA
##    for( i in 2:ns ) bweight[[i]] <- bw[1:i,i]
    if( nns > 1 ) for( i in 2:nns ) bweight[[i]] <- bw[1:i,i]

    raic <- array(z[[4L]], dim=c(ns,ns))
    aic <- list()
    aic[[1]] <- NA
##    for( i in 2:ns ) aic[[i]] <- raic[1:i,i]
    if( nns > 1 ) for( i in 2:nns ) aic[[i]] <- raic[1:i,i]

    a <- array(z[[5L]], dim=c(d,d,morder,ns))
    arcoef <- list()
##    for( i in 1:ns ) arcoef[[i]] <- array(a[,,,i], dim=c(d,d,morder))
    for( i in 1:nns ) arcoef[[i]] <- array(a[,,,i], dim=c(d,d,morder))

    e <- array(z[[6L]], dim=c(d,d,ns))
    v <- list()
##    for( i in 1:ns ) v[[i]] <- array(e[,,i], dim=c(d,d))
    for( i in 1:nns ) v[[i]] <- array(e[,,i], dim=c(d,d))

    aicbay <- z[[7L]][1:nns]
    start <- z[[8L]][1:nns]
    end <- z[[9]][1:nns]

##    blomar.out <- list( mean=z[[1L]], var=z[[2L]], bweight=bweight, aic=aic, arcoef=arcoef, v=v, eaic=z[[7L]],
##			init=z[[8L]], end=z[[9]] )
    blomar.out <- list( mean=z[[1L]], var=z[[2L]], bweight=bweight, aic=aic, arcoef=arcoef, v=v, eaic=aicbay,
			init=start, end=end )

    class( blomar.out ) <- "blomar"
    if( icflag == -1 ) cat(sprintf("\n ***** Warning : max.order is corrected n/(2*d) = %d\n\n", max.order))

    return( blomar.out )
}

print.blomar <- function(x, ...)
{
  id <- length(x$mean)
  cat("\n\n Mean ")
  for( i in 1:id ) cat(sprintf("  %f", x$mean[i]))
  cat("\n Variance ")
  for( i in 1:id ) cat(sprintf("  %f", x$var[i]))
  cat("\n")

  ns <- length(x$bweight)
  for( i in 1:ns ) {
    if( i != 1 ) {
      cat("\n\nAR-model fitted to\t! Bayesian weights\t! AIC with respect to the present data\n")
      cat("-------------------------------------------------------------------------------------\n")
      cat(sprintf("Current block\t\t! %f\t\t! %f\n", x$bweight[[i]][1], x$aic[[i]][1]))
      for( k in 1:(i-1) ) cat(sprintf("%i Period former block\t! %f\t\t! %f\n", k,x$bweight[[i]][k+1],x$aic[[i]][k+1]))
    }

    mf <- dim(x$arcoef[[i]])[3]
    cat("\n\n..........  Current model  ..........\n\n")
    cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n", x$init[i],x$end[i]))
    cat(sprintf(" aic = %f\n", x$eaic[i]))
    cat("\n Innovation variance matrix\n")
    print(x$v[[i]])
    cat(sprintf("\n AR coefficient matrix ( order %i )\n", mf))
    print(x$arcoef[[i]])
  }
}


bsubst <-
function (y, mtype, lag=NULL, nreg=NULL, reg=NULL, term.lag=NULL, cstep=5, plot=TRUE)
{
    if( mtype<0 || mtype>4 ) 
      stop("Allowed model type are
      1 : autoregressive model,
      2 : polynomial type non-linear model (lag's read in),
      3 : polynomial type non-linear model (lag's automatically set),
      4 : AR-model with polynomial mean value function\n")

# specification of number of regressors (mtype=1,2,4)
    if(is.null(nreg))  if( mtype != 3 )  stop("Number of regressors is not specified")

# specification of maximum time lag  (mtype=2)
    if(is.null(reg))  if( mtype == 2 )  stop("Regressors are not specified")

# specification of regressor (mtype=3)
    if(is.null(term.lag))  if( mtype == 3 ) stop("Maximum time lags are not specified")

    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum time lag use in the model
# number of regressors
    if( mtype == 3) {
      k <- 0
      lag5 <- term.lag[5]
      for( i in 1:lag5 )  k <- k + i*(i+1)/2 - 1
      nreg <- term.lag[1] + term.lag[2] + term.lag[3] + term.lag[4] + k
    }

    if(is.null(reg))  reg <- array(0, dim=c(3,nreg)) 
    if(is.null(term.lag))  term.lag <- rep(0,5)
    f <- ""                                          # specification of regressor (mtype=5)
    cnst <- 0                                        # constant value (mtype=6)

    z <- .Call("bsubst",
	as.double(y),
	as.integer(n),
	as.integer(mtype),
	as.integer(lag),
	as.integer(nreg),
	as.integer(cstep),
	as.integer(reg),
	as.integer(term.lag),
	as.character(f),
	as.double(cnst))

    mo <- z[[3L]]
    arcoefm=z[[6L]][1:mo]
    if( mtype != 4 ) {
      nps <- lag+1
      perr <- array(z[[23L]], dim=c(n,cstep))
      perr <- perr[nps:n,]
      lagh <- 100
      nn <- n-nps-1
      if( (lagh > nn) || (lagh == nn) ) lagh <- nn-1
      lag1 <- lagh + 1
      autcor1=z[[28L]][1:lag1]
    }


    if( plot == TRUE ) {

      if( mtype != 4 )  {
         mm <- 2
         if( mtype==1 ) mm <- mm+1
         nc <- (mm+cstep+1)/2
         par(mfcol=c(nc,2))
        for( i in 1:cstep ) {
          sig <- sqrt(z[[25L]][i])*0.5
#          hist( perr[(lag+1):n,i], breaks=c(-sig*10:1,0,sig*1:10), main="", xlab=paste(i,"-step ahead prediction error") )
          hist( perr[,i], breaks=c(-sig*10:1,0,sig*1:10), main="", xlab=paste(i,"-step ahead prediction error") )
        }
      }

      if( mtype == 1 ) {
        k <- lag
        k1 <- k+1
        plot((0:k), z[[9L]][1:k1], ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(m)-aicmin" )
        abline(h=0, lty=1)
      } else {
        plot((0:nreg), z[[9L]], ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(m)-aicmin" )
        abline(h=0, lty=1)
      }

      if( mtype==1 ) {
        x <- rep(0,121)
        for( i in 1:121 ) x[i] <- (i-1)/240
        plot(x, z[[29L]], type="l", xlab="Frequency", ylab=paste("Power Spectrum"))
      }

      if( mtype != 4 ) {
        plot( order<-c(0:lagh),autcor1,ylim=c(-1,1),bty="l",type="h", main="Autocorrelation of\n1-step ahead prediction error", xlab="Lines show +/-2sd\n ( sd = sqrt(1/n) )", ylab="peautcor" )
        abline( h=0, lty=1 )
        abline( h=2*sqrt(1/(n-lag)), lty=3 )
        abline( h=-2*sqrt(1/(n-lag)), lty=3 )
      }
      par(mfrow=c(1,1))
    }


    if(mtype == 1) {
          k <- lag
          k1 <- k+1
	bsubst.out <- list( ymean=z[[1L]], yvar=z[[2L]], v=z[[7L]][1:k1], aic=z[[8L]][1:k1], aicmin=z[[4L]], daic=z[[9L]][1:k1], order.maice=mo, v.maice=z[[5L]], arcoef.maice=arcoefm, v.bay=z[[11L]], aic.bay=z[[10L]], np.bay=z[[12L]], arcoef.bay=z[[13L]][1:k], ind.c=z[[14L]][1:k], parcor2=z[[15L]][1:k], damp=z[[16L]][1:k1], bweight=z[[17L]][1:k], parcor.bay=z[[18L]][1:k], eicmin=z[[19L]], esum=z[[20L]][1:k1], npmean=z[[21L]], npmean.nreg=z[[22L]], perr=perr, mean=z[[24L]], var=z[[25L]], skew=z[[26L]], peak=z[[27L]], peautcor=autcor1, pspec=z[[29]] ) }

    if(mtype==2 || mtype==3)
	bsubst.out <- list( ymean=z[[1L]], yvar=z[[2L]], v=z[[7L]], aic=z[[8L]], aicmin=z[[4L]], daic=z[[9L]], order.maice=mo, v.maice=z[[5L]], arcoef.maice=arcoefm, v.bay=z[[11L]], aic.bay=z[[10L]], np.bay=z[[12L]], arcoef.bay=z[[13L]], ind.c=z[[14L]], parcor2=z[[15L]], damp=z[[16L]], bweight=z[[17L]], parcor.bay=z[[18L]], eicmin=z[[19L]], esum=z[[20L]], npmean=z[[21L]], npmean.nreg=z[[22L]], perr=perr, mean=z[[24L]], var=z[[25L]], skew=z[[26L]], peak=z[[27L]], peautcor=autcor1 )

    if(mtype == 4)
	bsubst.out <- list( ymean=z[[1L]], yvar=z[[2L]], v=z[[7L]], aic=z[[8L]], aicmin=z[[4L]], daic=z[[9L]], order.maice=mo, v.maice=z[[5L]], arcoef.maice=arcoefm, v.bay=z[[11L]], aic.bay=z[[10L]], np.bay=z[[12L]], arcoef.bay=z[[13L]], ind.c=z[[14L]], parcor2=z[[15L]], damp=z[[16L]], bweight=z[[17L]], parcor.bay=z[[18L]], eicmin=z[[19L]], esum=z[[20L]], npmean=z[[21L]], npmean.nreg=z[[22L]] )

    return( bsubst.out )
}

exsar <-
function (y, max.order=NULL, plot=FALSE)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    morder <- max.order

    z <- .Call("exsar",
	as.double(y),
	as.integer(n),
	as.integer(morder))

    if( z[[12L]] != 0 ) stop("in FUNCT : SD is less than or equal to 0")

    if( plot == TRUE ) {
      plot((0:morder), z[[5L]], ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(m)-aicmin (Truncated at 40.0)")
      abline(h=0, lty=1) }

    m <- z[[6L]]
    exsar.out  <- list( mean=z[[1L]], var=z[[2L]], v=z[[3L]], aic=z[[4L]], aicmin=z[[7L]], daic=z[[5L]], order.maice=m, v.maice=z[[8L]], arcoef.maice=z[[9L]][1:m], v.mle=z[[10L]], arcoef.mle=z[[11L]][1:m] )
    return( exsar.out )
}


mlocar <-
function (y, max.order=NULL, span, const=0, plot=TRUE)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    morder <- max.order

    if( span < 1 ) span <- n
    ns <- as.integer((n-morder+span-1)/span)

    z <- .Call("mlocar",
	as.double(y),
	as.integer(n),
	as.integer(morder),
	as.integer(span),
	as.integer(const),
	as.integer(ns))

    a <- array(z[[3L]], dim=c(morder+const,ns))
    arcoef <- list()
    for(i in 1:ns) arcoef[[i]] <- a[1:z[[4L]][i],i]
    pspec <- array(z[[8L]], dim=c(121,ns))
    npre <- z[[9L]]
    order.const=z[[14L]]
    v.const=z[[15L]]
    aic.const=z[[16L]]
    npre[1] <- NA
    order.const[1] <- NA
    v.const[1] <- NA
    aic.const[1] <- NA

    if( plot == TRUE ) {
      oldpar <- par(no.readonly=TRUE)
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      ymin <- min(pspec)
      ymax <- max(pspec)
      par(mfrow=c(ns,1))
      nk <- ns
      if( ns == 4 ) {
        par(mfrow=c(2,2))
        nk <- 4
      }
      if( ns > 4 ) {
        par(mfrow=c(3,2))
        nk <- 6
      }
      if( ns > 6 ) {
        par(mfrow=c(3,3))
        nk <- 9
      }
      nn <- 0
      for(i in 1:ns) {
        if( nn == nk ) {
          par(ask=TRUE)
          nn <- 0
        }
        plot(x, pspec[,i], type="l", ylim=c(ymin,ymax), main=paste("y(", z[[6L]][i], "),...,y(", z[[7L]][i], ")"),
                          xlab="Frequency", ylab="Power Spectrum")
        nn <- nn+1
      }
      par(oldpar)
    }

    mlocar.out <- list( mean=z[[1L]], var=z[[2L]], ns=ns, order=z[[4L]], arcoef=arcoef, v=z[[5L]], init=z[[6L]], end=z[[7L]], 
			pspec=pspec, npre=npre, nnew=z[[10L]], order.mov=z[[11L]], v.mov=z[[12L]], aic.mov=z[[13L]],
			order.const=order.const, v.const=v.const, aic.const=aic.const )
    class( mlocar.out ) <- "mlocar"
    return( mlocar.out )
}

print.mlocar <- function(x, ...)
{

  cat(sprintf("\n\n Mean\t%f\n", x$mean))
  cat(sprintf(" Variance\t%f\n", x$var))

  ns <- x$ns
  for( i in 1:ns ) {
    if( i == 1 )
      cat(sprintf("\n Initial local model: (nnew = %i)\tvariance = %e\taic = %f\n", x$nnew[i], x$v.mov[i], x$aic.mov[i]))

    if( i != 1 ) {
      cat("\n\n >>>  The following two models are compared  <<<\n")
      np <- x$npre[i]+x$nnew[i]
      cat(sprintf(" Moving model: (npre = %i, nnew = %i)\tvariance = %e\taic = %f\n", x$npre[i], x$nnew[i], x$v.mov[i], x$aic.mov[i]))
      cat(sprintf(" Constant model: (npre+nnew = %i)\tvariance = %e\taic = %f\n", np, x$v.const[i], x$aic.const[i]))
      if( x$aic.mov[i] < x$aic.const[i] ) {
        cat(" --->  New model adopted\n")
      } else {
        cat(" --->  Constant model adopted\n")
      }
    }
    cat("\n\n..........  Current model  ..........\n\n")
    cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n", x$init[i],x$end[i]))
    cat(sprintf(" Innovation variance = %e\n\n", x$v[i]))
    cat(sprintf(" AR coefficients ( order %i ) \n", x$order[i]))
    print(x$arcoef[[i]])
  }

}


mlomar <-
function (y, max.order=NULL, span, const=0)
{
    n <- nrow(y)
    d <- ncol(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    icflag <- 0
    if( max.order > n/(2*d) )  icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))
    morder <- max.order

    calb <- rep(1,d)   # calibration for channel j (j=1,d)
    if( span < 1 ) span <- n
    ns <- as.integer((n-morder+span-1)/span)

    z <- .Call("mlomar",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.double(calb),
	as.integer(morder),
	as.integer(span),
	as.integer(const),
	as.integer(ns))

    nns <- z[[15L]]
    npre <- z[[3L]][1:nns]
    nnew <- z[[4L]][1:nns]
    order.mov <- z[[5L]][1:nns]
    aic.mov <- z[[6L]][1:nns]
    order.const <- z[[7L]][1:nns]
    aic.const <- z[[8L]][1:nns]
    order <- z[[9L]][1:nns]
    aic <- z[[10L]][1:nns]
    start <- z[[13L]][1:nns]
    end <- z[[14L]][1:nns]
    npre[1] <- NA
    order.mov[1] <- NA
    aic.mov[1] <- NA
    order.const[1] <- NA
    aic.const[1] <- NA

    a <- array(z[[11L]], dim=c(d,d,morder,ns))
    arcoef <- list()
##    for( i in 1:ns ) arcoef[[i]] <- array(a[,,,i],dim=c(d,d,z[[9L]][i]))
    for( i in 1:nns ) arcoef[[i]] <- array(a[,,,i],dim=c(d,d,z[[9L]][i]))

    e <- array(z[[12L]], dim=c(d,d,ns))
    v <- list()
##    for( i in 1:ns ) v[[i]] <- array(e[,,i],dim=c(d,d))
    for( i in 1:nns ) v[[i]] <- array(e[,,i],dim=c(d,d))

##    mlomar.out <- list( mean=z[[1L]], var=z[[2L]], ns=ns, order=z[[9L]], aic=z[[10L]], arcoef=arcoef, v=v,
##			init=z[[13L]], end=z[[14L]], npre=npre, nnew=nnew, order.mov=order.mov, aic.mov=aic.mov,
    mlomar.out <- list( mean=z[[1L]], var=z[[2L]], ns=nns, order=order, aic=aic, arcoef=arcoef, v=v,
			init=start, end=end, npre=npre, nnew=nnew, order.mov=order.mov, aic.mov=aic.mov,
			order.const=order.const, aic.const=aic.const )

    class( mlomar.out ) <- "mlomar"
    if( icflag == -1 ) cat(sprintf("\n ***** Warning : max.order is corrected n/(2*d) = %d\n\n", max.order))

    return( mlomar.out )
}

print.mlomar <- function(x, ...)
{

  id <- length(x$mean)
  cat("\n\n Mean ")
  for( i in 1:id ) cat(sprintf("  %f", x$mean[i]))
  cat("\n Variance ")
  for( i in 1:id ) cat(sprintf("  %f", x$var[i]))

  ns <- x$ns
  for( i in 1:ns ) {
    if( i == 1 )
      cat(sprintf("\n\n Initial local model: (nnew = %i), aic = %f\n", x$nnew[i],x$aic[i]))

    if( i != 1 ) {
      cat("\n\n >>>  The following two models are compared  <<<\n")
      np <- x$npre[i]+x$nnew[i]
      cat(sprintf(" Moving model: (npre = %i, nnew = %i)\taic = %f\n", x$npre[i], x$nnew[i], x$aic.mov[i]))
      cat(sprintf(" Constant model: (npre+nnew = %i)\taic = %f\n", np,x$aic.const[i]))
      if( x$aic.mov[i] < x$aic.const[i] ) {
        cat(" ---> New model adopted \n")
      } else {
        cat(" ---> Constant model adopted \n")
      }
    }
    cat("\n\n..........  Current model  ..........\n\n")
    cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n", x$init[i],x$end[i]))
    cat(sprintf(" aic = %f\n", x$aic[i]))
    cat("\n Innovation variance matrix\n")
    print(x$v[[i]])
    cat(sprintf("\n AR coefficient matrix ( order %i )\n", x$order[i]))
    print(x$arcoef[[i]])
  }

}


mulbar <-
function (y, max.order=NULL, plot=FALSE)
{
    n <- nrow(y)
    d <- ncol(y)
    calb<-rep(1,d)   # calibration of channel i (i=1,d)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of the order of AR-model
    icflag <- 0
    if( max.order > n/(2*d) )  icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))
    morder <- max.order

    z <- .Call("mulbar",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.double(calb),
	as.integer(morder))

    if( plot == TRUE ) {
      plot( (0:morder), z[[5L]], ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(m)-aicmin (Truncated at 40.0)" )
      abline(h=0, lty=1) }

    mulbar.out <- list( mean=z[[1L]], var=z[[2L]], v=z[[3L]], aic=z[[4L]], aicmin=z[[7L]], daic=z[[5L]],
	order.maice=z[[6L]], v.maice=z[[8L]], bweight=z[[9L]], integra.bweight=z[[10L]],
	arcoef.for=array(z[[11L]], dim=c(d,d,morder)), arcoef.back=array(z[[12L]], dim=c(d,d,morder)),
	pacoef.for=array(z[[13L]], dim=c(d,d,morder)), pacoef.back=array(z[[14L]], dim=c(d,d,morder)),
	v.bay=array(z[[15L]], dim=c(d,d)), aic.bay=z[[16L]] )

    if( icflag == -1 ) cat(sprintf("\n ***** Warning : max.order is corrected n/(2*d) = %d\n\n", max.order))

    return( mulbar.out )
}


mulmar <-
function (y, max.order=NULL, plot=FALSE) 
{
    n <- nrow(y)
    d <- ncol(y)
    calb <- rep(1, d)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    icflag <- 0
    if( max.order > n/(2*d) )  icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))

    z <- .Call("mulmar", as.double(y), as.integer(n), as.integer(d), 
        as.double(calb), as.integer(max.order))

    lag1 <- max.order+1
    v <- list()
    aic <- list()
    daic <- list()
    for (i in 1:d) {
        j <- ((i - 1) * lag1 + 1):(i * lag1)
        v[[i]] <- z[[3L]][j]
        aic[[i]] <- z[[4L]][j]
        daic[[i]] <- z[[5L]][j]
    }
    jnd <- list()
    subregcoef <- list()
    ind <- array(z[[10L]], dim = c(lag1 * d, d))
    a <- array(z[[11L]], dim = c(lag1 * d, d))
    for (i in 1:d) jnd[[i]] <- ind[(1:z[[9L]][[i]]), i]
    for (i in 1:d) subregcoef[[i]] <- a[(1:z[[9L]][[i]]), i]

    if (plot == TRUE) {
        par(mfrow=c(d,1))  
        for (i in 1:d) {
            plot((0:max.order), daic[[i]], ylim = c(0, 40), type = "l", 
                main = paste(" d=", i), xlab = "Lag", ylab = "AIC(m)-aicmin (Truncated at 40.0)")
            abline(h = 0, lty = 1)
        }
        par(mfrow=c(1,1))
    }

    mulmar.out <- list(mean = z[[1L]], var = z[[2L]], v = v, aic = aic, aicmin = z[[7L]],
        daic = daic, order.maice = z[[6L]], v.maice = z[[8L]], np = z[[9L]], jnd = jnd,
        subregcoef = subregcoef, rvar = z[[12L]], aicf = z[[13L]], respns = array(z[[14L]], dim = c(d, d)),
        matv = array(z[[16L]], dim=c(d,d)), morder = z[[18L]], arcoef = array(z[[17L]], dim=c(d,d,z[[18L]])), aicsum = z[[19L]])

    if( icflag == -1 ) cat(sprintf("\n ***** Warning : max.order is corrected n/(2*d) = %d\n\n", max.order))

    return(mulmar.out)
}


perars <-
function (y, ni, lag=NULL, ksw=0)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(ni))
    lag1 <- lag+1

    z <- .Call("perars",
	as.double(y),
	as.integer(n),
	as.integer(ni),
	as.integer(lag),
	as.integer(ksw))

    npara <- z[[3L]]
    ind <- array(z[[4L]], dim=c(lag1*ni+ksw,ni))
    jnd <- list()
    for( i in 1:ni )  jnd[[i]] <- ind[1:npara[i],i]

    a <- array(z[[5L]], dim=c(lag1*ni+ksw,ni))
    regcoef <- list()
    for( i in 1:ni )  regcoef[[i]] <- a[1:npara[i],i]

    perars.out <- list( mean=z[[1L]], var=z[[2L]], subset=jnd, regcoef=regcoef, rvar=z[[10L]], np=npara, 
                           aic=z[[6L]], v=array(z[[8L]], dim=c(ni,ni)), arcoef=array(z[[7L]], dim=c(ni,ni,z$morder)),
                           const=z[[9L]], morder=z[[11L]] )

    class( perars.out ) <- "perars"
    return( perars.out )

}

print.perars <- function(x, ...)
{
  cat(sprintf("\n\n Mean = %f\n", x$mean))
  cat(sprintf(" Variance = %f\n", x$var))

  ni <- nrow(x$v)
  for( i in 1:ni ) {
    cat(sprintf("\n ........  Regression Model for the regressand i = %i  ........\n", i))
    cat("\n subset")
    sorder <- length(x$subset[[i]])
    for( j in 1:sorder )  cat(sprintf("\t\t%i", x$subset[[i]][j]))
    cat("\n regcoef")
    for( j in 1:sorder ) cat(sprintf("\t%f", x$regcoef[[i]][j]))
    cat(sprintf("\n\n residual variance (rvar) = %f\n", x$rvar[i]))
    cat(sprintf(" number of parameter (np) = %i\n", x$np[i]))
    cat(sprintf(" AIC = n*log(rvar) + 2*np = %f\n\n", x$aic[i]))
  }

  cat("\n\nmatrix of regression coefficients\n")
  print(x$v)
  cat("\n\nregression coefficients within the present period\n")
  print(x$arcoef)
  cat("\n\nconstants within the regression models\n")
  print(x$const)
}


unibar <-
function (y, ar.order=NULL, plot=TRUE)
{
    n <- length(y)
    if( is.null(ar.order) ) ar.order <- as.integer(2*sqrt(n))
    ar.order1 <- ar.order+1

    z <- .Call("unibar",
	as.double(y),
	as.integer(n),
	as.integer(ar.order))

    if( plot == TRUE ) {
      oldpar <- par(no.readonly=TRUE)
      par(mfrow=c(3,1))
      plot((0:ar.order), z[[5L]], ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(m)-aicmin (Truncated at 40.0)")
      abline(h=0, lty=1)
      plot(z[[9L]], type="h", xlab="Lag", ylab="Partial autocorrelation")
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n),lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3)
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      plot(x, z[[17L]], type="l", xlab="Frequency", ylab="Power Spectral Density")
      par(oldpar)
      unibar.out <- list( mean=z[[1L]], var=z[[2L]], v=z[[3L]], aic=z[[4L]], aicmin=z[[7L]], order.maice=z[[6L]], v.maice=z[[8L]], bweight=z[[10L]][2:(ar.order1)], integra.bweight=z[[11L]], v.bay=z[[14L]], aic.bay=z[[13L]], np=z[[15L]], pacoef.bay=z[[12L]], arcoef=z[[16L]])
    } else {
      unibar.out <- list( mean=z[[1L]], var=z[[2L]], v=z[[3L]], aic=z[[4L]], aicmin=z[[7L]], daic=z[[5L]], 
order.maice=z[[6L]], v.maice=z[[8L]], pacoef=z[[9L]], bweight=z[[10L]][2:(ar.order1)], integra.bweight=z[[11L]],
v.bay=z[[14L]], aic.bay=z[[13L]], np=z[[15L]], pacoef.bay=z[[12L]], arcoef=z[[16L]],  pspec=z[[17L]] )
    }
    return( unibar.out )
}


unimar <-
function (y, max.order=NULL, plot=FALSE)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of AR-order
    morder <- max.order

    z <- .Call("unimar",
	as.double(y),
	as.integer(n),
	as.integer(morder)
)

    m <- z[[6L]]
    if( plot == TRUE ) {
      plot((0:morder), z[[5L]], ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(m)-aicmin (Truncated at 40.0)")
      abline(h=0, lty=1)
      unimar.out <- list( mean=z[[1L]], var=z[[2L]], v=z[[3L]], aic=z[[4L]], aicmin=z[[7L]], order.maice=m, v.maice=z[[8L]], arcoef=z[[9L]][1:m] )
    } else { 
    unimar.out <- list( mean=z[[1L]], var=z[[2L]], v=z[[3L]], aic=z[[4L]], aicmin=z[[7L]], daic=z[[5L]],  order.maice=m, v.maice=z[[8L]], arcoef=z[[9L]][1:m] )
    }
    return( unimar.out )
}


xsarma <-
function (y, arcoefi, macoefi)
{
    n <- length(y)
    arcoefi <- -arcoefi   # Initial estimates of AR coefficients
    macoefi <- -macoefi   # Initial estimates of MA coefficients
    p <- length(arcoefi)  # AR-ORDER
    q <- length(macoefi)  # MA-ORDER
    p01 <- c(arcoefi,macoefi)    # INITIAL ESTIMATES OF AR- AND  MA-COEFFICIENTS

    z <- .Call("xsarma",
	as.double(y),
	as.integer(n),
	as.integer(p),
	as.integer(q),
	as.double(p01))

    coef <- -z[[3L]]

    xsarma.out <- list( gradi=z[[1L]], lkhoodi=z[[2L]], arcoef=coef[1:p], macoef=coef[(p+1):(p+q)],
			grad=z[[4L]], alph.ar=z[[5L]], alph.ma=z[[6L]], lkhood=z[[7L]], wnoise.var=z[[8L]])
    return( xsarma.out )
}


#####   TIMSAC84   #####

decomp <- function(y, trend.order=2, ar.order=2, frequency=12, seasonal.order=1, log=FALSE, trade=FALSE,
	    diff=1, year=1980, month=1, miss=0, omax=99999.9, plot=TRUE)
{
    m1 <- trend.order
    m2 <- ar.order
    ilog <- 0
    if( log == TRUE ) ilog <- 1
    itrade <- 0
    if( trade == TRUE ) itrade <- 1

    n <- length(y)
    ipar <- rep(0, 9)
    ipar[1] <- trend.order
    ipar[2] <- ar.order
    ipar[3] <- frequency
    ipar[4] <- seasonal.order
    ipar[5] <- ilog
    ipar[6] <- itrade
    ipar[7] <- diff
    ipar[8] <- year
    ipar[9] <- month

    z <- .Call("decomp",
             as.double(y),
             as.integer(n),
             as.integer(ipar),
             as.integer(miss),
             as.double(omax))	

    ier <- z[[7L]]
    if( ier != 0 ) stop("Log-transformation cannot be applied to zeros and nagative numbers")

    aic=z[[6L]][1]
    lkhd=z[[6L]][2]
    sigma2=z[[6L]][3]
    tau1=z[[6L]][4]
    tau2=z[[6L]][5]
    tau3=z[[6L]][6]
    arcoef=z[[6L]][7:(6+m2)]
    tdf=z[[6L]][(7+m2):(13+m2)]

    if( plot == TRUE ) {
      oldpar <- par(no.readonly=TRUE)
      nc <- 2
      if( ar.order != 0 ) nc <- nc+1
      if( seasonal.order != 0 ) nc <- nc+1
      if( itrade == 1 ) nc <- nc+1
      if( nc > 3 ) par(mfrow=c((nc+1)/2,2))
      if( nc <= 3 ) par(mfrow=c(nc,1))

      if( ilog == 0 ) {
        matD <- array(c(y,z[[1L]]),dim=c(n,2))
        matplot(matD, pch=0, type="l", col=1:2, main="Original and Trend", xlab="", ylab="")
      } else {
        yy <- y
        for( i in 1:n )
          if( yy[i] > 0 ) {
            yy[i] <- log(yy[i])
          } else {
            yy[i] <- -1
          }
        ymax <- max(z[[1L]], yy) 
        ymin <- min(z[[1L]], yy)
        plot(0,0, type='n',xlim=c(0,n), ylim=c(ymin,ymax), xlab="", ylab="")
        y1 <- yy[1]
        for( i in 1:(n-1) ) {
          par(new=TRUE)
          y2 <- yy[i+1]
          if( (y1 != -1) && (y2 != -1) ) lines(c(i,i+1), c(y1,y2))
          y1 <- y2
        }
        par(new=TRUE)
        plot(z[[1L]], type='l', col=2, xlim=c(0,n), ylim=c(ymin,ymax), main="log(y) and Trend", xlab="", ylab="")
      }

      ymax <- max(z[[2L]], z[[3L]], z[[4L]], z[[5L]])
      ymin <- min(z[[2L]], z[[3L]], z[[4L]], z[[5L]])
      my <- max(ymax, abs(ymin))*1.5
      if( seasonal.order != 0 ) plot(z[[2L]], type="l", main= "Seasonal", xlab="", ylab="", ylim=c(-my,my))
      plot(z[[5L]], type="l", main= "Noise", xlab="", ylab="", ylim=c(-my,my))
      if( ar.order != 0 ) plot(z[[3L]], type="l", main="AR component", xlab="", ylab="", ylim=c(-my,my))
      if( itrade == 1) plot(z[[4L]], type="l", main="Trading Day Effect", xlab="", ylab="", ylim=c(-my,my))
      par(oldpar)
      decomp.out <- list(aic=aic, lkhd=lkhd, sigma2=sigma2, tau1=tau1, tau2=tau2, tau3=tau3, arcoef=arcoef, tdf=tdf)
    } else {
      decomp.out <- list(trend=z[[1L]], seasonal=z[[2L]], ar=z[[3L]], trad=z[[4L]], noise=z[[5L]],
                       aic=aic, lkhd=lkhd, sigma2=sigma2, tau1=tau1, tau2=tau2, tau3=tau3, arcoef=arcoef, tdf=tdf)
    }
    return(decomp.out)
}


baysea <- function(y, period=12, span=4, shift=1, forecast=0, trend.order=2, seasonal.order=1, year=0, month=1, out=0, rigid=1, zersum=1, delta=7, alpha=0.01, beta=0.01, gamma=0.1, spec=TRUE, plot=TRUE, separate.graphics=FALSE)
{
    if( seasonal.order > span ) stop("seasonal.order is smaller than or equal to span")
    if( span < 1 ) stop("span is greater than or equal to 1")
    if( trend.order < 1 ) stop("trend.order is greater than or equal to 1")

    ndata <- length(y)
    npf <- ndata+forecast

    ipara <- rep(0, 12)
    ipara[1] <- period
    ipara[2] <- span
    ipara[3] <- shift
    ipara[4] <- trend.order
    ipara[5] <- seasonal.order
    ipara[6] <- 0           # logt
    ipara[7] <- year
    ipara[8] <- month
    ipara[9] <- 1                            # nday
    if( spec == TRUE )  ipara[10] <- 1       # spectrum estimation option
    if( spec == FALSE ) ipara[10] <- 0
    ipara[11] <- out                        # ioutd : outlier correction option 

    para <- rep(0,8)
    para[1] <- rigid     # controls the rigidity of the seasonal component
    para[2] <- 1         # wtrd
    para[3] <- 1         # dd
    para[4] <- zersum    # controls the sum of the seasonals within a period
    para[5] <- delta     # controls the leap year effect
    para[6] <- alpha     # controls prior variance of initial trend
    para[7] <- beta      # controls prior variance of initial seasonal
    para[8] <- gamma     # controls prior variance of initial sum of seasonal

    arft <- rep(0, 3)
    arfs <- rep(0, 3)
    arfn <- rep(0, 3)

#   subroutine arcoef ---> subroutine partar
    iart <- 0
#    for(i in 1:3) if(arft[i] != 0.0) iord <- i
#    if(iord != 0) {
#      arm <- matrix(0, dim=c(3,3))
#      for( i in 1:3 ) arm[i,i] <- arft[i]
#      arm[2,1] <- arm[1,1]-arft[2]*arm[1,1]
#      arm[3,1] <- arm[2,1]-arft[3]*arm[2,2]
#      arm[3,2] <- arm[2,2]-arft[3]*arm[2,1]
#      for( i in 1:iord ) arft[i] <- arm[iord,i]
#    }
    iars <- 0
    iarn <- 0

    is <- period*seasonal.order
    iprd <- 2
    if( period == 1 ) iprd <- 1
    lftrn <- trend.order + iart
    lfsea <- (seasonal.order+iars)*period + iarn
    idc <- lftrn*iprd + 1
    idcx <- lfsea*2 + 1
    if( period>1 && idc<idcx ) idc <- idcx
    if( period>1 && idc<period*2-1 ) idc <- period*2-1
    ipara[12] <- idc

    z <- .Call("baysea",
             as.double(y),
             as.integer(ndata),
             as.integer(forecast),
             as.integer(ipara),
             as.double(para),
             as.double(arft),
             as.double(arfs),
             as.double(arfn),
             as.integer(iart),
             as.integer(iars),
             as.integer(iarn))	

    outlier <- NULL
    if( out != 0 ) outlier <- z[[1L]]
    tre <- z[[3L]]
    sea <- z[[4L]]
    tday <- NULL
    if( year != 0 ) tday <- z[[5L]]
    irr <- z[[6L]]
    adj <- z[[7L]]
    est <- z[[8L]]
    psds <- z[[9L]]
    psdt <- z[[10L]]
    abic <- z[[11L]]

    baysea.out <- list(outlier=outlier, trend=tre, season=sea, tday=tday, irregular=irr, adjust=adj, smoothed=est, aveABIC=abic)

    if( spec == TRUE ) {
      z1 <- spec.baysea(y, period, tre, trend.order, sea, seasonal.order, irr, adj)
      spec1 <- z1$irregular.spec
      spec2 <- z1$adjusted.spec
      spec3 <- z1$differenced.trend
      spec4 <- z1$differenced.season
      if( plot == TRUE ) plot.baysea(y, period, outlier, tre, trend.order, sea, seasonal.order, tday, irr, adj, est, psdt, psds, spec1, spec2, spec3, spec4, spec, separate.graphics )
      ir.spec <- list(acov=spec1$acov, acor=spec1$acor, mean=spec1$mean, v=spec1$v, aic=spec1$aic, parcor=spec1$parcor, rspec=spec1$rspec) 
      ad.spec  <- list(acov=spec2$acov, acor=spec2$acor, mean=spec2$mean, v=spec2$v, aic=spec2$aic, parcor=spec2$parcor, rspec=spec2$rspec) 
      diff.trend <- list(acov=spec3$acov, acor=spec3$acor, mean=spec3$mean, v=spec3$v, aic=spec3$aic, parcor=spec3$parcor)
      diff.season <- list(acov=spec4$acov, acor=spec4$acor, mean=spec4$mean, v=spec4$v, aic=spec4$aic, parcor=spec4$parcor)
      baysea.out <- c(baysea.out, irregular.spec=ir.spec, adjusted.spec=ad.spec, differenced.trend=diff.trend, differenced.season=diff.season)
    } else if ( plot == TRUE ) {
      ir.spec <- NULL
      ad.spec <- NULL
      diff.trend <- NULL
      diff.season <- NULL
      plot.baysea(y, period, outlier, tre, trend.order, sea, seasonal.order, tday, irr, adj, est, psdt, psds, ir.spec, ad.spec, diff.trend, diff.season, spec, separate.graphics)
    }

    return( baysea.out )
}

spec.baysea <- function(y, period, trend, trend.order, season, seasonal.order, irregular, adjust)
{
    ndata <- length(y)
    npf <- length(trend)

    lag <- min((ndata-1), 60)
    lag1 <- lag+1
#    ifpl <- min(3*sqrt(ndata), 50, lag)
    ifpl <- min(30, ndata-1)
    ifpl1 <- ifpl+1

 # SPECTRUM OF IRREGULAR                               
    mode <- 1
    z <- .Call("spgrh",
             as.double(irregular),
             as.integer(ndata),
             as.integer(lag1),
             as.integer(ifpl1),
             as.integer(mode),
             as.integer(period))

    pxx <- z[[7L]]
    sxx <- rep(0, lag1)
    for( i in 1:lag1 ) {
      t <- abs(pxx[i])
      sxx[i] <- log10(t)
    }

    if( z[[8L]] == 2600 ) cat("\n ***** Warning : Accuracy of computation lost\n\n" )

    irregular.spec <- list(n=ndata, acov=z[[1L]], acor=z[[2L]], mean=z[[3L]], v=z[[4L]], aic=z[[5L]], parcor=z[[6L]], rspec=pxx, rpspec=sxx)

# SPECTRUM OF DIFFERENCED ADJUSTED SERIES
    n1 <- ndata-1
    dadj <- adjust
    for( i in 1:n1 ) dadj[i] <- dadj[i+1]-dadj[i]
    dadj <- dadj[1:n1]
    lag <- min((n1-1), 60)
    lag1 <- lag+1

    z <- .Call("spgrh",
             as.double(dadj),
             as.integer(n1),
             as.integer(lag1),
             as.integer(ifpl1),
             as.integer(mode),
             as.integer(period))

    pxx <- z[[7L]]
    sxx <- rep(0, lag1)
    for( i in 1:lag1 ) {
      t <- abs(pxx[i])
      sxx[i] <- log10(t)
    }

    if( z[[8L]] == 2600 ) cat(" ***** Warnning : Accuracy of computation lost\n" )

    adjusted.spec <- list(n=n1, acov=z[[1L]], acor=z[[2L]], mean=z[[3L]], v=z[[4L]], aic=z[[5L]], parcor=z[[6L]], rspec=pxx, rpspec=sxx)

#  PARCOR OF TREND.ORDER TIME(S) DIFFERENCED TREND SERIES
    n1 <- ndata		
    trendd <- trend
    for( j in 1:trend.order ) {
      n1 <- n1-1
      for( i in 1:n1 ) trendd[i] <- trendd[i+1] - trendd[i]
    }
    trendd <- trendd[1:n1]
    lag <- min((n1-1), 60)
    lag1 <- lag+1
    mode <- 0

    z <- .Call("spgrh",
             as.double(trendd),
             as.integer(n1),
             as.integer(lag1),
             as.integer(ifpl1),
             as.integer(mode),
             as.integer(period))

    if( z[[8L]] == 2600 ) cat("\n ***** Warning : Accuracy of computation lost\n\n" )

    differenced.trend <- list(n=n1, acov=z[[1L]], acor=z[[2L]], mean=z[[3L]], v=z[[4L]], aic=z[[5L]], parcor=z[[6L]])

#  PARCOR OF SEASONAL.ORDER TIME(S) DIFFERENCED SEASONAL SERIES
    n1 <- ndata
    seasond <- season
    for( j in 1:seasonal.order ) {
      n1 <- n1-period
      for( i in 1:n1 ) seasond[i] <- seasond[i+period] - seasond[i]
    }
    seasond <- seasond[1:n1]

    lag <- min((n1-1), 60)
    lag1 <- lag+1

    z <- .Call("spgrh",
             as.double(seasond),
             as.integer(n1),
             as.integer(lag1),
             as.integer(ifpl1),
             as.integer(mode),
             as.integer(period))

    if( z[[8L]] == 2600 ) cat("\n ***** Warning : Accuracy of computation lost\n\n" )

    differenced.season <- list(n=n1, acov=z[[1L]], acor=z[[2L]], mean=z[[3L]], v=z[[4L]], aic=z[[5L]], parcor=z[[6L]])

    spec.baysea.out <- list(irregular.spec=irregular.spec, adjusted.spec=adjusted.spec, differenced.trend=differenced.trend, differenced.season=differenced.season)
}


plot.baysea <- function(y, period, outlier, trend, trend.order, season, seasonal.order, tday, irregular, adjust, smoothed, psdt, psds,
                        irregular.spec, adjusted.spec, differenced.trend, differenced.season, spec, separate.graphics)

{
    oldpar <- par(no.readonly=TRUE)
    ndata <- length(y)
    npf <- length(trend)

    ymax1 <- max(y, trend, adjust, smoothed)
    ymin1 <- min(y, trend, adjust, smoothed)

    plot(y, type="l", main="Original Data", xlab="time", xlim=c(0,npf), ylim=c(ymin1,ymax1))
    nw <- 1

    if( separate.graphics == TRUE ) {
      dev.new()
    } else {
      par(ask=TRUE)
    }

    plot(trend, type="l", main="Trend and 2*(post SD)", cex.main=0.9,  xlab="time", ylab="", xlim=c(0,npf), ylim=c(ymin1,ymax1))
    par(new=TRUE)
    xtem <- trend + psdt
    plot(xtem, type="l", lty=3, main="", xlab="", ylab="", xlim=c(0,npf), ylim=c(ymin1,ymax1))
    par(new=TRUE)
    xtem <- trend - psdt
    plot(xtem, type="l", lty=3, main="", xlab="", ylab="", xlim=c(0,npf), ylim=c(ymin1,ymax1))

    if( separate.graphics == TRUE )  dev.new()
    plot(adjust, pch=0, type="l", main="Adjusted = Original Data - Seasonal - Trading.Day.Comp - Outlier", cex.main=0.9,  xlab="time", ylab="", xlim=c(0,npf), ylim=c(ymin1,ymax1))

    if( separate.graphics == TRUE )  dev.new()
    plot(smoothed, pch=0, type="l", main="Smoothed = Trend + Seasonal + Trading.Day.Comp", cex.main=0.9,  xlab="time", ylab="", xlim=c(0,npf), ylim=c(ymin1,ymax1))

    ymax2 <- max(irregular)
    ymin2 <- min(irregular)
    if( seasonal.order != 0 ) {
      ymax2 <- max(season, ymax2)
      ymin2 <- min(season, ymin2)
    }
    if( is.null(tday) == FALSE) {
      ymax2 <- max(tday, ymax2)
      ymin2 <- min(tday, ymin2)
    }
    my <- max(ymax2, abs(ymin2))*1.5
    if( seasonal.order != 0 ) {
      if( separate.graphics == TRUE )  dev.new()
      plot(season, type="l", main= "Seasonal and 2*(post SD)", cex.main=0.9,  xlab="time", ylab="", ylim=c(-my,my))
      par(new=TRUE)
      xtem <- season + psds
      plot(xtem, type="l", lty=3, main="", xlab="", ylab="", ylim=c(-my,my))
      par(new=TRUE)
      xtem <- season - psds
      plot(xtem, type="l", lty=3, main="", xlab="", ylab="", ylim=c(-my,my))
    }

    if( separate.graphics == TRUE )  dev.new()
    par(mfrow=c(2,1))
    plot(irregular, type="l", main="Irregular = Original Data - Trend - Seasonal - Trading.Day.Comp", cex.main=0.9,  xlab="time", ylab="", ylim=c(-my,my))
    vy <- sd(irregular)*5
    plot(irregular, type="l", main="Irregular ( Scaled by the Standard Deviation)", cex.main=0.9,  xlab="time", ylab="", ylim=c(-vy,vy))
    par(mfrow=c(1,1))

 # SPECTRUM OF IRREGULAR                               
    if( spec == TRUE ) {
      lag <- length(irregular.spec$acov)-1
      n <- irregular.spec$n
      if( separate.graphics == TRUE )  dev.new()
      par(mfrow=c(3,1))
      plot((0:lag), irregular.spec$acor, type="h", main="Autocorrelation & Parcor of Irregular",  ylab="Autocorrelation", xlab="lag", ylim=c(-1,1))
      plot(irregular.spec$parcor, type="h", ylab="Parcor", xlab="lag", ylim=c(-1,1))
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n), lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3)

      it <- irregular.spec$rpspec * 10
      ymin <- (min(irregular.spec$rpspec) - 1) * 10
      ymax <- (max(irregular.spec$rpspec) + 1) * 10
      plot(it, type='l', main="High Order AR-Spectrum as an approximation to periodgram ( Order is fixed at 30 )", ylab="Rational Spectrum", xlab="", ylim=c(ymin,ymax))
      par(mfrow=c(1,1))

# SPECTRUM OF DIFFERENCED ADJUSTED SERIES
      lag <- length(adjusted.spec$acov)-1
      n <- adjusted.spec$n
      if( separate.graphics == TRUE )  dev.new()
      par(mfrow=c(3,1))
      plot((0:lag), adjusted.spec$acor, type="h", main="Autocorrelation & Parcor of Differenced Adjusted Series", ylab="Autocorrelation", xlab="lag", ylim=c(-1,1))
      plot(adjusted.spec$parcor, type="h", ylab="Parcor", xlab="lag", ylim=c(-1,1))
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n), lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3)

      it <- adjusted.spec$rpspec * 10
      ymin <- (min(adjusted.spec$rpspec) - 1) * 10
      ymax <- (max(adjusted.spec$rpspec) + 1) * 10
      plot(it, type='l', main="High Order AR-Spectrum as an approximation to periodgram( Order is fixed at 30 )", xlab="", ylab="Rational Spectrum", ylim=c(ymin,ymax))
      par(mfrow=c(1,1))

#  PARCOR OF TREND.ORDER TIME(S) DIFFERENCED TREND SERIES
      lag <- length(differenced.trend$acov)-1
      n <- differenced.trend$n
      if( separate.graphics == TRUE )  dev.new()
      par(mfrow=c(2,1))
      plot((0:lag), differenced.trend$acor, type="h", main=paste( trend.order, "Time(s) Differenced Trend Series" ), ylab="Autocorrelation", cex.main=0.9,  xlab="lag", ylim=c(-1,1))
      plot(differenced.trend$parcor, type="h", ylab="Parcor", xlab="lag", ylim=c(-1,1))
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n), lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3)
      par(mfrow=c(1,1))

#  PARCOR OF SEASONAL.ORDER TIME(S) DIFFERENCED SEASONAL SERIES
     lag <- length(differenced.season$acov)-1
      n <- differenced.season$n
      if( separate.graphics == TRUE )  dev.new()
      par(mfrow=c(2,1))
      plot((0:lag), differenced.season$acor, type="h", main=paste( seasonal.order, "time(s) Differenced Seasonal Series"), cex.main=0.9,  ylab="Autocorrelation", xlab="lag", ylim=c(-1,1))
      plot(differenced.season$parcor, type="h", ylab="Parcor", xlab="lag", ylim=c(-1,1))
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n), lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3)
    }

    par(oldpar)
}

.noGenerics <- TRUE

options(warn.FPU=FALSE)
