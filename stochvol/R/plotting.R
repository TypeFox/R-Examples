paradensplot <- function(x, showobs = TRUE, showprior = TRUE, showxlab = TRUE,
		       	 mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
 if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
 if (!is.logical(showobs)) stop("If provided, argument 'showobs' must be TRUE or FALSE.")
 if (!is.null(simobj)) {
  if (!is(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
  sim <- TRUE
 } else sim <- FALSE
 oldpar <- par(mar=mar)
 paranames <- c(quote(mu), quote(phi), quote(sigma), quote(nu))
 cutat1 <- c(FALSE, TRUE, FALSE, FALSE)
 for (i in 1:ncol(x$para)) {
  mydensplot(x$para[,i], show.obs=showobs, main=paste("Density of", paranames[i]),
	     cutat1=cutat1[i], showxlab=showxlab, mgp = mgp, ...)
  if (isTRUE(showprior)) {
   paras <- x$priors[[i]]
   vals <- seq(from=par('usr')[1], to=par('usr')[2], len=1000)
   if (i == 1) lines(vals, dnorm(vals, paras[1], paras[2]), col=8, lty=2)
   else if (i == 2) lines(vals, .5*dbeta((vals+1)/2, paras[1], paras[2]), col=8, lty=2)
   else if (i == 3) lines(vals, 2*dnorm(vals, 0, sqrt(paras[1])), col=8, lty=2)
   else if (i == 4) lines(vals, dunif(vals, x$priors$nu[1], x$priors$nu[2]), col = 8, lty = 2)
   if (sim && (i <= 3 || length(simobj$para) == 4)) {
    points(simobj$para[i], 0, col = 3, cex = 2, pch = 16)
   }
  }
 }
 par(oldpar)
 invisible(x)
}

paratraceplot <- function(x, mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
 if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
 if (!is.null(simobj)) {
  if (!is(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
  sim <- TRUE
 } else sim <- FALSE
 oldpar <- par(mar=mar)
 paranames <- c(quote(mu), quote(phi), quote(sigma), quote(nu))
 for (i in 1:ncol(x$para)) {
  mytraceplot(x$para[,i], xlab="", mgp = mgp,
	      main=paste("Trace of ", paranames[i], " (thinning = ", x$thinning$para,")", sep=''), ...)
  if (sim && (i <= 3 || length(simobj$para) == 4)) {
   abline(h = simobj$para[i], col = 3, lty = 2)
  }
 }
 par(oldpar)
 invisible(x)
}

volplot <- function(x, forecast = 0, dates = NULL, show0 = FALSE,
		    col = NULL, forecastlty = NULL, tcl = -.4,
		    mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), simobj = NULL, ...) {
 if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
 if (!is.null(simobj)) {
  if (!is(simobj, "svsim")) stop("If provided, simobj must be an 'svsim' object.")
  sim <- TRUE
 } else sim <- FALSE
 oldpar <- par(mar = mar)
 where <- grep("%", dimnames(x$summary$latent)[[2]])
 obj <- t(100*exp(x$summary$latent[,where,drop=FALSE]/2))  # monotone transformation!
 qs <- dim(obj)[1]
 timelen <- dim(obj)[2]
 if (is.null(qs) | all(is.na(obj))) stop("No quantiles to plot.")
 if (is.null(col)) {
  cols <- rep(8, qs)
  cols[dimnames(obj)[[1]] == "50%"] <- 1
 } else cols <- col
 if (is.null(forecastlty)) forecastlty <- 2
 
 if (is(forecast, "svpredict") | (is.numeric(forecast) & length(forecast) == 1 & all(forecast != 0))) { # also draw future values
  thintime <- x$thinning$time
  
  if (thintime != 1) {
   lasth <- as.integer(gsub("h_", "", dimnames(x$latent)[[2]][dim(x$latent)[2]]))
   if (length(x$y) > lasth) {
       warning(paste("Thinning for time 'thintime' has not been set to one during sampling. This means we are forecasting conditional on h_", lasth, " and not on h_", length(x$y), ".", sep=''))
   }
  }
  
  if(is.numeric(forecast) & length(forecast) == 1 & all(forecast >= 1)) {
#   warning("Calling prediction method.")
   forecast <- predict(x, forecast)
  }
  if(!is(forecast, "svpredict")) stop("Argument 'forecast' must be a single nonnegative integer, or of class type 'svpredict'.")
  
  futlen <- dim(forecast)[2]
  
  xs <- matrix(rep(seq(timelen, timelen + futlen/thintime, len=futlen+1), qs), nrow=futlen+1)
  quants <- as.numeric(gsub("%", "", dimnames(obj)[[1]]))/100
  ys <- rbind(obj[,timelen], t(matrix(apply(100*exp(forecast/2), 2, quantile, quants), nrow=qs)))
 
  if (futlen/thintime > .01*timelen) {  # increase xlim to give space for forecast
   if (thintime == 1) {
    xlim <- c(0, timelen+futlen/thintime)
   } else {
    xlim <- c(1, timelen+futlen/thintime)
   }
  } else {
   xlim <- NULL
  }
 } else xlim <- NULL
 
 if(exists("sd", x$summary)) {
  mymain <- paste("Estimated conditional volatilities in percent (", paste(dimnames(obj)[[1]], collapse=' / '),
		 " posterior quantiles)", sep = '')
 } else {
  mymain <- paste("Estimated volatilities in percent (", paste(dimnames(obj)[[1]], collapse=' / '),
		 " posterior quantiles)", sep = '')
 }

 ts.plot(t(obj), gpars=list(xlim=xlim, col=cols, xlab='', xaxt='n', mgp=mgp, tcl=tcl,
			    main = mymain, ...))

 if (sim) {
  lines(100*simobj$vol, col = 3)
 }
 
 if (is(forecast, "svpredict")) {
  for (i in 1:qs) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
 }
 
 ax <- axis(1, tick=FALSE, labels=FALSE)  # just automagic axis ticks, don't draw yet

 if (show0) { # also draw latent0:
  thintime <- x$thin$time
  xs <- matrix(rep(c(1-1/thintime,1), qs), nrow=2)
  where <- grep("%", names(x$summary$latent0))
  ys <- rbind(100*exp(x$summary$latent0[where]/2), obj[,1])
  for (i in 1:qs) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
 }
   
 if (is.null(dates)) {
  dates <- c(0L, as.integer(gsub("h_", "", dimnames(x$latent)[[2]])))
  if (max(ax) > length(dates)) {  # means we are probably forecasting and need extra axis labels
   dates <- c(dates, seq(length(dates), max(ax), by=dates[2]-dates[1]))
  }
 } else {
  if (is(dates, "Date")) dates <- as.character(dates)
  if (length(dates) != ncol(x$latent)) {
   stop("Length of argument 'dates' differs from ncol(x$latent).")
  }
  dates <- c('', dates)
  ax <- ax[ax != 0]  # avoid "zero" tick
 }
 axis(1, at=ax, labels=dates[ax+1], mgp=mgp, tcl=tcl)

 if(exists("sd", x$summary)) {
  ts.plot(100*x$summary$sd[,where], gpars=list(xlim=xlim, col=cols, xlab='', xaxt='n', mgp=mgp, tcl=tcl,
			        main = paste("Estimated volatilities in percent (",
				       paste(dimnames(obj)[[1]], collapse=' / '),
				       " posterior quantiles)", sep=''), ...))
  
  if (sim) {
   standardizer <- sqrt(simobj$para$nu / (simobj$para$nu - 2))
   lines(100*simobj$vol*standardizer, col = 3)
  }
  
  if (is(forecast, "svpredict")) {
   standardizer <- sqrt(x$para[,"nu"] / (x$para[,"nu"] - 2))
   ys <- rbind(100*x$summary$sd[timelen,where,drop=FALSE],
               t(matrix(apply(100*exp(forecast/2)*standardizer, 2, quantile, quants), nrow=qs)))

   for (i in 1:qs) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
  }
  axis(1, at=ax, labels=dates[ax+1], mgp=mgp, tcl=tcl)
 }

 par(oldpar)
 invisible(x)
}

plot.svdraws <- function(x, forecast = NULL, dates = NULL,
			 show0 = FALSE, showobs = TRUE, showprior = TRUE, col = NULL,
			 forecastlty = NULL, tcl = -0.4,
			 mar = c(1.9, 1.9, 1.7, .5), mgp = c(2, .6, 0),
			 simobj = NULL, ...) {
 oldpar <- par(mfrow=c(1,1))
 if (ncol(x$para) == 4) {
  layout(matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10), 4, byrow = TRUE))
 } else {
  layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), 3, byrow = TRUE))
 }
 volplot(x, dates = dates, show0 = show0, forecast = forecast,
	 forecastlty = forecastlty, col = col, tcl = tcl, mar = mar,
	 mgp = mgp, simobj = simobj, ...)
 paratraceplot(x, mar = mar, mgp = mgp, simobj = simobj, ...)
 paradensplot(x, showobs = showobs, showprior = showprior,
	      showxlab = FALSE, mar = mar, mgp = mgp, simobj = simobj, ...)
 par(oldpar)
 invisible(x)
}

# modified density plot (from coda package)
mydensplot <- function(x, show.obs = TRUE, bwf, main = "", ylim, cutat1=FALSE, showxlab=TRUE, mgp = c(2,.6,0), tcl=-.4, ...) 
{
    xx <- as.matrix(x)
    for (i in 1:nvar(x)) {
        y <- xx[, i, drop = TRUE]
        if (missing(bwf)) 
            bwf <- function(x) {
                x <- x[!is.na(as.vector(x))]
                return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
            }
        bw <- bwf(y)
        width <- 4 * bw
        if (max(abs(y - floor(y))) == 0 || bw == 0) 
            hist(y, prob = TRUE, main = main, ...)
        else {
            scale <- "open"
        if (isTRUE(cutat1)) {
	     if (1-max(y) < 2*bw) {
	      scale <- "cutat1"
	      y <- c(y, 2 - y)
	      if (1+min(y) < 2*bw) {
	       scale <- "cutatboth"
	       y <- c(y, -2 - y, 2 - y)
	      }
	     } else if (1+min(y) < 2*bw) {
	      scale <- "cutat-1"
	      y <- c(y, -2 - y)
	     }
	    }
	    else if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
            	if (min(y) >= 0 && min(y) < 2 * bw) {
                  scale <- "proportion"
                  y <- c(y, -y, 2 - y)
                }
            }
            else if (min(y) >= 0 && min(y) < 2 * bw) {
                scale <- "positive"
                y <- c(y, -y)
            }
	    else scale <- "open"
            dens <- density(y, width = width)
            if (scale == "proportion") {
                dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 
                  1]
                dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
            }
            else if (scale == "positive") {
                dens$y <- 2 * dens$y[dens$x >= 0]
                dens$x <- dens$x[dens$x >= 0]
            }
	    else if (scale == "cutat1") {
	    	dens$y <- 2 * dens$y[dens$x <= 1]
	    	dens$x <- dens$x[dens$x <= 1]
	    }
	    else if (scale == "cutat-1") {
	    	dens$y <- 2 * dens$y[dens$x >= -1]
	    	dens$x <- dens$x[dens$x >= -1]
	    }
	    else if (scale == "cutatboth") {
	    	dens$y <- 3 * dens$y[dens$x >= -1 & dens$x <= 1]
	    	dens$x <- dens$x[dens$x >= -1 & dens$x <= 1]
	    }
	    if (missing(ylim)) 
                ylim <- c(0, max(dens$y))
            
	    plot(dens, ylab = "", main = main, type = "l", 
		  ylim = ylim, xlab="", mgp = mgp, tcl = tcl, ...)
            if(isTRUE(showxlab)) {
	       if (is.R()) {
                  mtext(paste("N =", niter(x), "  Bandwidth =",
			      formatC(dens$bw)), side=1, line=2.7, cex=.7)
               } else {
                  mtext(paste("N =", niter(x), "  Bandwidth =",
			      formatC(bw)), side=1, line=2.7, cex=.7)
               }
	    }
            if (show.obs) 
                lines(y[1:niter(x)], rep(max(dens$y)/100, niter(x)), 
                  type = "h")
        }
        if (!is.null(varnames(x)) && is.null(list(...)$main)) 
            title(paste("Density of", varnames(x)[i]))
    }
    return(invisible(x))
}

# modified traceplot (from coda)
mytraceplot <- function (x, smooth = FALSE, col = 1:6, type = "l", ylab = "", xlab = "Iterations", mgp = c(2,.6,0), tcl = -.4, ...) 
{
    x <- mcmc.list(x)
    args <- list(...)
    for (j in 1:nvar(x)) {
        xp <- as.vector(time(x))
        yp <- if (nvar(x) > 1) 
            x[, j, drop = TRUE]
        else x
        yp <- do.call("cbind", yp)
        matplot(xp, yp, xlab = xlab, ylab = ylab, type = type, 
            col = col, mgp = mgp, tcl = tcl, ...)
        if (!is.null(varnames(x)) && is.null(list(...)$main)) 
            title(paste("Trace of ", varnames(x)[j], " (thin = ", attr(x, "thinning")$thinpara,")", sep=''))
        if (smooth) {
            scol <- rep(col, length = nchain(x))
            for (k in 1:nchain(x)) lines(lowess(xp, yp[, k]), 
                col = scol[k])
        }
    }
}

plot.svresid <- function(x, origdata = NA,
			 mains = c("Residual plot", "Q-Q plot"),
			 mar = c(2.9, 2.7, 2.2, .5),
			 mgp = c(1.7, .6, 0), ...) {
 
 if (any(is.na(origdata))) {
  oldpar <- par(mfrow=c(1, 2), mar=mar, mgp=mgp)
 } else {
  oldpar <- par(mfrow=c(2, 2), mar=mar, mgp=mgp)
  plot.default(origdata, ylab='Original values', xlab='Time', xaxt='n', ylim=c(-1,1)*max(abs(origdata)), main="Original data", ...)
  where <- seq(1, length(origdata), length=min(7, length(origdata)))
  axis(1, at = where, labels = names(origdata)[where])
  qqnorm(origdata, main=paste(mains[2], "for original data"))
  qqline(origdata, probs = c(0.01, 0.99))
 }

 # NEW: Cater for conditional t-distributions
 if (!is.null(attr(x, "nu"))) {
  terr <- TRUE
  nu <- attr(x, "nu")
  xlab <- paste("Theoretical quantiles from a t-distribution with", round(nu, 2), "df")
 } else {
  terr <- FALSE
  nu <- Inf
  xlab <- "Theoretical quantiles from a standard normal distribution"
 }
 
 plot.default(x, ylab = paste("M", substring(attr(x, "type"), 2), ' standardized residuals', sep = ""),
	      xlab='Time', xaxt='n', ylim=c(-1,1)*max(abs(x)),
	      main=mains[1], ...)
 
 if (!terr) abline(h=qnorm(c(.025, .975)), lty=2)
 where <- seq(1, length(x), length=min(7, length(x)))
 axis(1, at = where, labels = gsub("r_", "", names(x)[where]))
 qqplot(qt(ppoints(length(x)), df = nu), x,
	main=paste(mains[2], "for", attr(x, "type"), "standardized residuals"),
	xlab = xlab, ylab = "Sample quantiles")
 qqline(x, probs = c(0.01, 0.99), distribution = function(x) qt(x, df = nu))
 par(oldpar)
 invisible(x)
}
