lines.muhaz <- function(x, ...) {
    lines(x$est.grid, x$haz.est, ...)
}

muhaz <- function(times, delta, subset, min.time, max.time, bw.grid, bw.pilot,
                  bw.smooth, bw.method="local", b.cor="both", n.min.grid=51,
                  n.est.grid=101, kern="epanechnikov")
{
    method <- pmatch(bw.method, c("global", "local", "knn"))
    if (is.na(method))
        stop("\nbw.method MUST be one of: 'global', 'local', or 'knn'\n")

    b.cor <- pmatch(b.cor, c("none", "left", "both"))
    if (is.na(b.cor))
        stop("\nb.cor MUST be one of: 'none', 'left', or 'both'\n")
    b.cor <- b.cor - 1

    kern <- pmatch(kern, c("rectangle", "epanechnikov", "biquadratic",
                           "triquadratic"))
    if (is.na(kern))
           stop("\n kern MUST be one of: 'rectangle', 'epanechnikov', 'biquadratic', or 'triquadratic'\n")
    kern <- kern - 1

    if ( missing(times) )
        stop("\nParameter times is missing\n")

    nobs <- length(times)

    if ( missing(delta) ) {
           delta <- rep(1, nobs)
       } else {
           if ( length(delta) != nobs ) {
               stop("\ntimes and delta MUST have the same length\n")
           }
       }

    if ( missing(subset) ) {
        subset <- rep(TRUE, nobs)
    } else {
        if ( is.logical(subset) ) {
            if ( length(subset) != nobs ) {
                     stop("\ntimes and subset MUST have the same length\n")
                 }
        } else {
            stop("\nsubset MUST contain ONLY logical values (T or F)\n")
        }
    }
      # pick only the observations in subset
    times <- times[subset]
    delta <- delta[subset]

      # sort data, on times
    ix <- order(times)
    times <- times[ix]
    delta <- delta[ix]

    nobs <- length(times)

    if ( missing(min.time) ) {
        startz <- 0
    } else {

        if (min.time > times[1]) {
            warning("minimum time > minimum Survival Time\n\n")
        }

        startz <- min.time
    }

    if ( missing(max.time) ) {
    # use the time corresponding to 10 survivors
        sfit <- survfit( Surv(times,delta) ~ 1 )
        endz <- approx(sfit$n.risk,sfit$time,xout=10)$y
    } else {
        if (max.time > times[nobs]) {
            warning("maximum time > maximum Survival Time\n\n")
            endz <- times[nobs]
        } else {
            endz <- max.time
        }
    }
    if (startz > endz)
        stop("\nmin.time MUST be < max.time\n")
    # use default values as proposed by Mueller (bw.pilot and bw.smooth)
    nz <- sum(delta)
    if ( missing(bw.pilot) ) {
        bw.pilot <- endz/8/(nz^.2)
    }

    if ( missing(bw.smooth) ) {
        bw.smooth <- 5 * bw.pilot
    }

      # build the minimization grid
    z <- seq(startz, endz, len=n.min.grid)

      # build the estimation grid
    zz <- seq(startz, endz, len=n.est.grid)

      # this are irrelevant now
    endl <- startz
    endr <- endz

      # whatever the method, these are common input parameters

    pin.common <- list(times=times, delta=delta,
                       nobs=nobs, min.time=startz, max.time=endz,
                       n.min.grid=n.min.grid, min.grid=z,
                       n.est.grid=n.est.grid,
                       bw.pilot=bw.pilot,  bw.smooth= bw.smooth,
                       method=method, b.cor=b.cor, kernel.type=kern)

    if (method < 3) {
        if (missing(bw.grid)) {
            bw.grid <- seq(0.2*bw.pilot, 20*bw.pilot, len=25)
        }
        gridb <- length(bw.grid)
        m1 <- method - 1
        ans <- .Fortran("newhad",
                        as.integer(nobs),
                        as.double(times),
                        as.integer(delta),
                        as.integer(kern),
                        as.integer(m1),
                        as.double(z),
                        as.integer(n.min.grid),
                        as.double(zz),
                        as.integer(n.est.grid),
                        as.double(bw.pilot),
                        as.double(bw.grid),
                        as.integer(gridb),
                        as.double(endl),
                        as.double(endr),
                        as.double(bw.smooth),
                        as.integer(b.cor),
                        fzz = double(n.est.grid),
                        bopt = double(n.min.grid),
                        bopt1 = double(n.est.grid),
                        msemin = double(n.min.grid),
                        biasmin = double(n.min.grid),
                        varmin = double(n.min.grid),
                        imsemin = double(1),
                        globlb = double(1),
                        globlmse = double(gridb),
			PACKAGE = "muhaz")

        if (method == 1) {
            ans <- list(pin=pin.common, est.grid=zz, haz.est=ans$fzz,
                        imse.opt=ans$imsemin, bw.glob=ans$globlb,
                        glob.imse=ans$globlmse, bw.grid=bw.grid)
        } else {
            ans <- list(pin=pin.common, est.grid=zz, haz.est=ans$fzz,
                        bw.loc=ans$bopt, bw.loc.sm=ans$bopt1, msemin=ans$msemin,
                        bias.min=ans$biasmin, var.min=ans$varmin,
                        imse.opt=ans$imsemin)
        }
    } else if (method == 3) {
        kmin <- 2
        kmax <- as.integer(nz/2)
        m1 <- 2
        ans <- .Fortran("knnhad",
                        as.integer(nobs),
                        as.double(times),
                        as.integer(delta),
                        as.integer(kern),
                        as.integer(m1),
                        as.integer(n.min.grid),
                        as.double(z),
                        as.integer(n.est.grid),
                        as.double(zz),
                        as.double(bw.pilot),
                        as.double(endl),
                        as.double(endr),
                        as.double(bw.smooth),
                        as.integer(b.cor),
                        fzz = double(n.est.grid),
                        kopt = as.integer(kmin),
                        as.integer(kmax),
                        bopt = double(n.min.grid),
                        bopt1 = double(n.est.grid),
                        kimse = double(kmax-kmin+1),
			PACKAGE = "muhaz")
        ans <- list(pin=pin.common, est.grid=zz, haz.est=ans$fzz,
                    bw.loc=ans$bopt, bw.loc.sm=ans$bopt1, k.grid=kmin:kmax,
                    k.imse= ans$kimse, imse.opt=ans$kimse[ans$kopt-kmin+1],
                    kopt=ans$kopt)
    } else {
        stop("\nmethod MUST be 1, 2, or 3\n")
    }
    class(ans) <- "muhaz"
    ans
}

pehaz <- function(times, delta = NA, width = NA, min.time = 0, max.time = NA)
{
    if(is.na(delta[1]))
        delta <- rep(1,length(times))

    if(is.na(max.time))
        max.time <- max(times)

    if(is.na(width)) {
        nu <- sum(delta)
        width <- (max.time - min.time)/(8 * nu^0.2)
    }
    cat("\nmax.time=", max.time)
    cat("\nwidth=", width)

    nbins <- ceiling((max.time - min.time)/width)
    cat("\nnbins=", nbins)
    cat("\n")
    cuts <- rep(NA, nbins + 1)
    hazfun <- rep(NA, nbins)
    fusum <- rep(NA, nbins)
    evnum <- rep(NA, nbins)
    atrisk <- rep(NA, nbins)
    cuts[1] <- min.time

    for(i in 1:nbins) {
        left <- min.time + (i - 1) * width
                                        # bins are considered left-continuous
        right <- min.time + i * width
        cuts[i + 1] <- right
        ind <- (times < right) & !(times < left)
	# ind indicates pts terminating within interval
        fu <- times[!(times < left)]
        fu[(fu > right)] <- right
        fu <- fu - left
        sumfu <- sum(fu)
        numevnt <- sum(delta[ind])
        haz <- numevnt/sumfu
        fusum[i] <- sumfu
        evnum[i] <- numevnt
        atrisk[i] <- sum(!(times < left))
        hazfun[i] <- haz
    }

                                        # save call
    obj.call <- match.call()

      	# organize results
        # evnum = number of events in each bin
        # atrisk = number at risk in each bin
        # fusum = sum of f/u time in each bin

    result <- list(call = obj.call, Width = width , Cuts = cuts,
                   Hazard = hazfun, Events = evnum, At.Risk = atrisk,
                   F.U.Time = fusum)

    class(result) <- "pehaz"
    result
}

plot.muhaz <- function(x, ylim, type, xlab, ylab, ...) {

      y <- x$haz.est
      if(missing(ylim))
	ylim<-c(0,max(y))
      if(missing(type))
        type<-'l'
      if(missing(xlab))
        xlab<-'Follow-up Time'
      if(missing(ylab))
        ylab<-'Hazard Rate'
      plot(x$est.grid, y, type,  ylim=ylim, xlab=xlab, ylab=ylab, ...)

      return (invisible())
}

plot.pehaz <- function(x, xlab = "Time", ylab = "Hazard Rate", ...)
{
    lenh <-length(x$Hazard)
    hvals<-c(x$Hazard, x$Hazard[lenh])
    plot(x$Cuts, hvals, type = "s", ylab = ylab,
         xlab = xlab, ...)
}

lines.pehaz <- function(x, lty = 2, ...)
{
  lenh <-length(x$Hazard)
  hvals<-c(x$Hazard, x$Hazard[lenh])
  lines(x$Cuts, hvals, type="s", lty = lty, ...)
}

print.pehaz <- function( x , ... )
  {
    cat("\nCall:\n")
    print(x$call)

    cat("\nBin Width:\n")
    print(x$Width)

    cat("\nCuts Defining the Bins:\n")
    print(x$Cuts)

    cat("\nHazard Estimate for Each Bin:\n")
    print(x$Hazard)

    cat("\nNumber of Events in Each Bin:\n")
    print(x$Events)

    cat("\nNumber at Risk in Each Bin:\n")
    print(x$At.Risk)

    cat("\nTotal Follow-up Time in Each Bin:\n")
	print(x$F.U.Time)

    invisible()
}
summary.muhaz <- function(object,...) {

  cat("\nNumber of Observations ..........", object$pin$nobs)
  cat("\nCensored Observations ...........", object$pin$nobs-sum(object$pin$delta))

  cat("\nMethod used ..................... ")
  y <- object$pin$method
  if (y == 1) {
    cat("Global")
  } else if (y == 2) {
    cat("Local")
  } else if (y == 3) {
    cat("Nearest Neighbor")
  } else {
           cat("Unknown method")
         }

  cat("\nBoundary Correction Type ........ ")

  y <- object$pin$b.cor
  if (y == 0) {
    cat("None")
  } else if (y == 1) {
    cat("Left Only")
  } else if (y == 2) {
    cat("Left and Right")
  } else {
    cat("Unknown")
      }

  cat("\nKernel type ..................... ")
  y <- object$pin$kernel.type
  if (y == 0) {
    cat("Rectangle")
  } else if (y == 1) {
    cat("Epanechnikov")
  } else if (y == 2) {
    cat("Biquadratic")
  } else if (y == 3) {
    cat("Triquadratic")
  } else {
    cat("Unknown")
  }

  cat("\nMinimum Time ....................", round(object$pin$min.time,2))
  cat("\nMaximum Time ....................", round(object$pin$max.time,2))

  cat("\nNumber of minimization points ...", object$pin$n.min.grid)
  cat("\nNumber of estimation points .....", object$pin$n.est.grid)

  cat("\nPilot Bandwidth .................", round(object$pin$bw.pilot,2))
  cat("\nSmoothing Bandwidth .............", round(object$pin$bw.smooth,2))

  y <- object$pin$method
  if (y == 1) {
    cat("\nOptimal Global Bandwidth ........", round(object$bw.glob,2))
  } else if (y == 3) {
    cat("\nOptimal Nearest Neighbor ........", object$kopt)
  }

  cat("\nMinimum IMSE ....................", round(object$imse.opt,2))
  cat("\n")

  return (invisible())
}
kphaz.fit <- function(time, status, strata, q = 1, method = "nelson")
{
########################################################################
#   "kphaz.fit" is an S function which calculates a K-M-type hazard
#   function estimate.
#   Required Arguments:
#   time: vector of time values; all values must be greater than
#         or equal to zero. Missing values (NAs) are allowed.
#   status: vector of status values.  The values are 0 for
#           censored or 1 for uncensored (dead). Missing values
#           (NAs) are allowed. Must have the same length as time.
#   Optional Arguments:
#   strata: an optional vector that will be used to divide the
#           subjects into disjoint groups. Each group generates a
#           hazard curv.  Missing values (NAs) are allowed.
#           If missing, all subjects are assumed to be in the
#           same strata.
#   q: The number of failure times combined. Default is 1.
#   method: Type of hazard estimation made.
#           "product-limit" - estimator of the cumulative hazard
#               at the ordered failure times, t[j]; j=1,...,J
#               H[t[j]] = sum(t[i] <= t[j]) -log(1 - status[i]/(n-i+1))
#               Also provides an estimate of the variance
#               Var(H[t[j]]) = sum(t[i] <= t[j]) status[i]/((n-i+1)*(n-i))
#           "nelson" - estimator of the cumulative hazard
#               at the ordered failure times, t[j]; j=1,...,J
#               H[t[j]] = sum(t[i] <= t[j]) status[i]/(n-i+1)
#               Also provides an estimate of the variance
#               Var(H[t[j]]) = sum(t[i] <= t[j]) status[i]/((n-i+1)^2)
#   Returns:
#   time: a vector containing the death times at which estimates
#         were made
#   haz: a vector containing the hazard estimate at each time
#        in time.
#   var: a vector containing the variance estimates of the hazard
#        estimate at each time. Only returned if method is "Nelson".
#   strata: a vector which divides the hazard estimate into
#           disjoint groups. This vector is returned only if
#           'strata' is defined when 'kphaz.fit' is called.
#   Method:
#   (1) Estimate H[t[j]] and Var(H[t[j]]) at each of the distinct
#       death times according to the method selected.
#       (a) For the "nelson" method:
#               H[t[j]] = sum(t[i] <= t[j]) status[i]/(n-i+1)
#               Var(H[t[j]]) = sum(t[i] <= t[j]) status[i]/((n-i+1)^2)
#       (b) For the "product-limit" method:
#               H[t[j]] = sum(t[i] <= t[j]) -log(1 - status[i]/(n-i+1))
#               Var(H[t[j]]) = sum(t[i] <= t[j]) status[i]/((n-i+1)*(n-i))
#   (2) Define the hazard estimate at time (t[q+j]+t[j])/2 to be
#       haz[(t[q+j]+t[j])/2] = (H[t[q+j]]-H[t[j]])/(t[q+j]-t[j])
#       and the variance
#       var[(t[q+j]+t[j])/2] =
#            (Var(H[t[q+j]])-Var(H[t[j]]))/(t[q+j]-t[j])^2
########################################################################
########################
# Check Input Arguments
########################
#
# time
#
	if(missing(time)) stop("Argument \"time\" is missing, with no default")
	if(any(is.na(time)))
		stop("Time values can not be NA")
	if(any(is.nan(time)))
		stop("Time values can not be Infinite")
	if(any(time < 0))
		stop("Time values must be >= 0")
	if(any(!is.numeric(time))) stop("Time must be a numeric vector")	#
# Status
#
	if(missing(status))
		stop("Argument \"status\" is missing, with no default")
	if(any(is.na(status)))
		stop("Status values can not be NA")
	if(any(!is.numeric(status)))
		stop("Status must be a numeric vector")
	if(length(status) != length(time))
		stop("No. of observations in \"time\" and \"status\" must match"
			)
	status <- as.integer(status)
	status[status != 0] <- 1
	if(all(status == 0)) stop("No events occur in this data set")	#
# strata
#
	if(missing(strata))
		qstrata <- FALSE
	else {
		if(length(strata) != length(time))
			stop("\"Strata\" vector is the wrong length")
		qstrata <- TRUE
	}
	if(!is.numeric(q))
		stop("Agument \"q\" must be a numberic value")
	if(is.na(q))
		stop("q may not be NA")
	if(is.nan(q))
		stop("q may not be Infinite")
	q <- as.integer(q)
	if(q < 1) stop("q must be positive")	#
# Method
#
	imethod <- pmatch(method, c("nelson", "product-limit"))
	if(is.na(imethod)) stop("method must be one of \"nelson\" or \"product-limit\""
			)	##########################
# Sort the data by strata
##########################
	if(!qstrata)
		strata <- rep(1, length(time))
	ind <- order(strata)
	strata <- strata[ind]
	time <- time[ind]
	status <- status[ind]
	ustrata <- unique(strata)
	##################################################
# Initalize the output variables
# dtime = death times at which estimates are made
# haz = hazard estimate at dtime
# var = variance of hazard estimate at dtime
# dstrata = strata indicator
##################################################
	dtime <- vector()
	haz <- vector()
	var <- vector()
	dstrata <- vector()	#############################################
# Find hazard estimates for each unique
# strata.
#############################################
	for(j in 1:length(ustrata)) {
#
# Define the values of time and status which are included in
# the current strata
#
		cur.strata <- ustrata[j]
		cur.time <- time[strata == cur.strata]
		cur.status <- status[strata == cur.strata]	#
# Sort the current values of time and status by time
#
		ind <- order(cur.time)
		cur.time <- cur.time[ind]
		cur.n <- length(cur.time)
		cur.status <- cur.status[ind]	#
# Define cur.dtime = the list of unique death times
#
		cur.dtime <- unique(cur.time[cur.status != 0])
		cur.nd <- length(cur.dtime)
		if(cur.nd > q) {
#
# Calculate the current stratas hazard estimates
# and variance estimates at each distinct hazard time.
#
			H <- vector()
			VH <- vector()
			for(i in 1:cur.nd) {
				temp.status <- cur.status[cur.time <= cur.dtime[
				  i]]
				temp.n <- 1:length(temp.status)
	# Calculate the hazard.
				if(imethod == 1)
				  H[i] <- sum(temp.status/(cur.n - temp.n + 1))
				else if(imethod == 2)
				  H[i] <- sum( - log(1 - (temp.status/(cur.n -
				    temp.n + 1))))
				if(imethod == 1)
				  VH[i] <- sum(temp.status/((cur.n - temp.n + 1
				    ) * (cur.n - temp.n + 1)))
				else VH[i] <- sum(temp.status/((cur.n - temp.n +
				    1) * (cur.n - temp.n)))
			}
#
# Calculate the new death times and estimates based on q.
#
			for(i in seq(1, cur.nd - q)) {
				dtime <- c(dtime, ((cur.dtime[q + i] +
				  cur.dtime[i])/2))
				ttt <- (cur.dtime[q + i] - cur.dtime[i])
				haz <- c(haz, (H[q + i] - H[i])/ttt)
				var <- c(var, (VH[q + i] - VH[i])/(ttt^2))
				dstrata <- c(dstrata, cur.strata)
			}
		}
	}
#########################################################
# Return time,haz,var and strata if required
#########################################################
	time <- dtime
	strata <- dstrata
	if(qstrata)
		return(list(time=time, haz=haz, var=var, strata=strata))
	else return(list(time=time, haz=haz, var=var))
}
kphaz.plot <- function(fit, ...)
{
########################################################################
#   "kphaz.plot" is an S function plots the K-M-type hazard
#   function estimate generated by kphaz.fit.
#   Required Arguments:
#   fit: results of a call to kphaz.fit.
#   ...: additional arguments for the plot function
########################################################################
#
# Make sure time and haz exist and are the same length
#
	if(any(names(fit) == "time") & (any(names(fit) == "haz"))) {
		time <- fit$time
		haz <- fit$haz
	}
	else stop("Argument \"fit\" must be the result of a call to \"kphaz.fit\""
			)	#
	if(length(time) != length(haz)) stop(
			"Argument \"fit\" must be the result of a call to \"kphaz.fit\""
			)	#
# Check to see if there are any strata
#
	qstrata <- any(names(fit) == "strata")
	if(qstrata)
		strata <- fit$strata
	else strata <- rep(1, length(time))
	if(length(strata) != length(haz)) stop(
			"Argument \"fit\" must be the result of a call to \"kphaz.fit\""
			)	#
# Define, ustrata, the number of unique strata
#
	ustrata <- unique(strata)
	good <- 1:length(ustrata)
	for(i in 1:length(ustrata)) {
		cur.strata <- ustrata[i]
		ind <- strata == ustrata[i]
		if(all(is.na(haz[ind] | is.nan(haz[ind]))))
			good <- good[good != i]
	}
	ustrata <- ustrata[good]	#
# Check to see if there are any plots to be made
#
	if(length(ustrata) < 1) stop("No plots")	#
# Make the first plot
#
	ind <- (strata == ustrata[1]) & (!is.nan(haz))
	xmax <- max(time)
	ymax <- max(haz[((!is.nan(haz)) & (!is.na(haz)))])
	x <- time[ind]
	y <- haz[ind]
	if(min(x) > 0) {
		x <- c(0, x)
		y <- c(0, y)
	}
	plot(x, y, xlim = c(0, xmax), ylim = c(0, ymax), xlab = "Time",
		ylab = "Hazard", type = "s", ...)	#
# Make any subsequent plots
#
	if(length(ustrata) > 1) {
		for(i in 2:length(ustrata)) {
			ind <- strata == ustrata[i]
			x <- time[ind]
			y <- haz[ind]
			if(min(x) > 0) {
				x <- c(0, x)
				y <- c(0, y)
			}
			lines(stepfun(x, y), lty = i)
		}
	}
#
# Return
#
	invisible()
}
