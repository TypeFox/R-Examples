"plot.potd" <- 
function(x, ...)
{
    rawdata <- x$data
    n <- length(as.numeric(rawdata))
    times <- attributes(rawdata)$times
    if(is.character(times) || inherits(times, "POSIXt") ||
       inherits(x, "date") || inherits(x, "dates")) {
        times <- as.POSIXlt(times)
        gaps <- as.numeric(difftime(times[2:n], times[1:(n-1)],
            units = "days")) * x$intensity
    }
    else gaps <- as.numeric(diff(times)) * x$intensity
    data <- as.numeric(rawdata)
    threshold <- x$threshold
    par.ests <- x$par.ests
    xi <- par.ests[1]
    beta <- par.ests[4]
    residuals <- logb(1 + (xi * (data - threshold))/beta)/xi
    choices <- c("Point Process of Exceedances", "Scatterplot of Gaps",
		 "Qplot of Gaps", "ACF of Gaps", "Scatterplot of Residuals",
		 "Qplot of Residuals", "ACF of Residuals", "Go to GPD Plots")
    tmenu <- paste("plot:", choices)
    pick <- 1
    lastcurve <- NULL
    while(pick > 0) {
        pick <- menu(tmenu, title = 
		     "\nMake a plot selection (or 0 to exit):")
        if(pick %in% 1:7) lastcurve <- NULL
	switch(pick,
            {
	       plot(times, rawdata, type = "h", sub = paste("Point process of",
	         length(as.numeric(rawdata)), "exceedances of threshold",
                 format(signif(threshold, 3))), ...)  
            },
	    {
	      plot(gaps, ylab = "Gaps", xlab = "Ordering", ...)
	      lines(lowess(1:length(gaps), gaps))
	    },
	      qplot(gaps, ...),
	      acf(gaps, lag.max = 20, ...),
	    {
	      plot(residuals, ylab = "Residuals", xlab = "Ordering", ...)
	      lines(lowess(1:length(residuals), residuals))
	    },
	      qplot(residuals, ...),
	      acf(residuals, lag.max = 20, ...),
	      lastcurve <- plot.gpd(x, ...))
    }
    invisible(lastcurve)
}

"pot" <- 
function(data, threshold = NA, nextremes = NA, run = NA,
    picture = TRUE, ...)
{
    n <- length(as.numeric(data))
    times <- attributes(data)$times
    if(is.null(times)) {
        times <- 1:n
        attributes(data)$times <- times
        start <- 1
        end <- n
        span <- end - start
    }
    else {
        start <- times[1]
        end <- times[n]
        span <- as.numeric(difftime(as.POSIXlt(times)[n],
            as.POSIXlt(times)[1], units = "days"))
    }
  
    if(is.na(nextremes) && is.na(threshold))
       	stop("Enter either a threshold or the number of upper extremes")
    if(!is.na(nextremes) && !is.na(threshold))
	stop("Enter EITHER a threshold or the number of upper extremes")
    if(!is.na(nextremes))
	threshold <- findthresh(as.numeric(data), nextremes)
    if(threshold > 10) {
	factor <- 10^(floor(log10(threshold)))
	cat(paste("If singularity problems occur divide data",
                  "by a factor, perhaps", factor, "\n"))
    }
    exceedances.its <- structure(data[data > threshold], times =
        times[data > threshold])
    n.exceed <- length(as.numeric(exceedances.its))
    p.less.thresh <- 1 - n.exceed/n
    if(!is.na(run)) {
       	exceedances.its <- decluster(exceedances.its, run, picture)
       	n.exceed <- length(exceedances.its)
    }
    intensity <- n.exceed/span
    exceedances <- as.numeric(exceedances.its)
    xbar <- mean(exceedances) - threshold
    s2 <- var(exceedances)
    shape0 <- -0.5 * (((xbar * xbar)/s2) - 1)
    extra <- ((length(exceedances)/span)^( - shape0) - 1)/shape0
    betahat <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
    scale0 <- betahat/(1 + shape0 * extra)
    loc0 <- 0
    theta <- c(shape0, scale0, loc0)
    negloglik <- function(theta, exceedances, threshold, span)
    {
        if((theta[2] <= 0) || (min(1 + (theta[1] * (exceedances -
            theta[3])) / theta[2]) <= 0))
	    f <- 1e+06
	else {
	    y <- logb(1 + (theta[1] * (exceedances - theta[3])) / theta[2])
	    term3 <- (1/theta[1] + 1) * sum(y)
	    term1 <- span * (1 + (theta[1] * (threshold - theta[3])) /
                         theta[2])^(-1/theta[1])
	    term2 <- length(y) * logb(theta[2])
	    f <- term1 + term2 + term3
	}
	f
    }
    fit <- optim(theta, negloglik, hessian = TRUE, ..., exceedances =
                 exceedances, threshold = threshold, span = span)
    if(fit$convergence)
        warning("optimization may not have succeeded")
    par.ests <- fit$par
    varcov <- solve(fit$hessian)
    par.ses <- sqrt(diag(varcov))   
    beta <- par.ests[2] + par.ests[1] * (threshold - par.ests[3])
    par.ests <- c(par.ests, beta)
    out <- list(n = length(data), period = c(start, end), data = 
        exceedances.its, span = span, threshold = threshold,
        p.less.thresh = p.less.thresh, n.exceed = n.exceed, run = run,
	par.ests = par.ests, par.ses = par.ses, varcov = varcov, 
	intensity = intensity, nllh.final = fit$value, converged
	= fit$convergence)
    names(out$par.ests) <- c("xi", "sigma", "mu", "beta")
    names(out$par.ses) <- c("xi", "sigma", "mu")
    class(out) <- "potd"
    out
}

