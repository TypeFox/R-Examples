"turnogram" <-
function(series, intervals=c(1, length(series)/5), step=1, complete=FALSE, two.tailed=TRUE, FUN=mean, plotit=TRUE, level=0.05, lhorz=TRUE, lvert=FALSE, xlog=TRUE) {
	call <- match.call()
	data <- deparse(substitute(series))
	fun <- deparse(substitute(FUN))
	if (is.null(class(series)) || class(series) != "ts")
		stop("series must be a single regular time series")
	Unit <- attr(series, "units")
	UnitTxt <- GetUnitText(series)
	# Determine if there are many successive sequences of zeros
	X <- as.vector(series)
	N <- length(X)
	strips <- c(X[1]-1, X[1:(N-1)]) != X | X != 0
  	StrippedSeries <- X[strips]
  	# Test the length of the series, range and step...
  	n <- length(X)
  	ns <- length(StrippedSeries)
  	if (n > ns + 0.1 * n)
  		warning("There are more than 10% of data as successive zeros, turnogram could be biased!")
    if (length(intervals) < 2)
    	stop("Interval must be a vector with 2 values: (min, max)")
    if (intervals[1] < 1 || intervals[1] > n/3)
    	stop("Interval must be larger or equal to 1, and smaller or equal to n/3")
    if (intervals[2] < 1 || intervals[2] > n/3)
    	stop("Interval must be larger or equal to 1, and smaller or equal to n/3")	
    Inter.vec <- seq(intervals[1], intervals[2], step)
    if (length(Inter.vec) < 2)
    	stop("Less than 2 intervals. Redefine intervals or step")
	n.vec <- nturns.vec <- nturns.min.vec <- nturns.max.vec <- Info.vec <- Info.min.vec <- Info.max.vec <- Inter.vec
    
    # Calculate the first step for the turnogram
	turnogram.step1 <- function(Series, Interval, Complete, Fun) {
		if (Complete == FALSE) {		# We just start intervals from the first observation
			if (length(Series) %% Interval !=0) Series <- c(Series, rep(NA, Interval - (length(Series) %% Interval)))
			dim(Series) <- c(Interval, length(Series) %/% Interval)
			x <- apply(Series, 2, Fun, na.rm=TRUE)		# Because there is much chance to get some NAs appended at the end of the series!
			n <- length(x)
			Nturns <- turnpoints(x)$nturns
			res <- list(n=n, nturns=Nturns)
		} else {					# We calculate each possible interval
			n <- Nturns.vec <- NULL
			for (j in Interval:1) {
				Ser <- Series[j:length(Series)]
				if (length(Ser) %% Interval !=0) Ser <- c(Ser, rep(NA, Interval - (length(Ser) %% Interval)))
				dim(Ser) <- c(Interval, length(Ser) %/% Interval)
				x <- apply(Ser, 2, Fun, na.rm=TRUE)
				Nturns.vec[j] <- turnpoints(x)$nturns	
				n[j] <- length(x)
			}
			res <- list(n=n, nturns=Nturns.vec)
		}
		res
	}
    
    # Calculate all first steps (n and nturns) for the turnogram
	if (complete == FALSE) {
		for (i in 1:length(Inter.vec)) {
			res <- turnogram.step1(X, Inter.vec[i], Complete=complete, Fun=FUN)
			n.vec[i] <- res$n
			nturns.vec[i] <- res$nturns
		}
		# Calculate I (bits of information) according to either Gleissberg (n <= 50), or normal approximations (n > 50)
		Info.vec <- -log(pgleissberg(n.vec, nturns.vec, two.tailed=two.tailed), base=2)
		if (two.tailed == TRUE) {	# We have to change sign if k > mu.
			rightpart <- nturns.vec > 2 * (n.vec - 2) / 3
			Info.vec[rightpart] <- -Info.vec[rightpart]
		}
		# By default the extraction level is set to the interval corresponding to the maximum info value
		Level <- Inter.vec[match(max(Info.vec), Info.vec)]
		res <- list(interval=Inter.vec, n=n.vec, turns=nturns.vec, info=Info.vec, level=Level)
	} else {
		for (i in 1:length(Inter.vec)) {
			res <- turnogram.step1(X, Inter.vec[i], Complete=complete, Fun=FUN)
			n.vec[i] <- max(res$n)		# To change this!!!
			nturns <- res$nturns
			nturns.vec[i] <- mean(nturns)
			nturns.min.vec[i] <- min(nturns)
			nturns.max.vec[i] <- max(nturns)
			infos <- -log(pgleissberg(res$n, nturns, two.tailed=two.tailed), base=2)
			if (two.tailed == TRUE) {	# We have to change sign if k > mu.
				rightpart <- nturns > 2 * (res$n - 2) / 3
				infos[rightpart] <- -infos[rightpart]
			}
			Info.vec[i] <- mean(infos)
			Info.min.vec[i] <- min(infos)
			Info.max.vec[i] <- max(infos)
			# By default the extraction level is set to the interval corresponding to the maximum info value
			Level <- Inter.vec[match(max(Info.vec), Info.vec)]
			res <- list(interval=Inter.vec, n=n.vec, turns=nturns.vec, turns.min=nturns.min.vec, turns.max=nturns.max.vec, info=Info.vec, info.min=Info.min.vec, info.max=Info.max.vec, level=Level)
		}
	}
	as.data.frame(res)
	res$call <- call
	res$data <- data
	if (complete == TRUE) res$type <- "Complete" else res$type <- "Simple"
	res$fun <- fun
	if (two.tailed == TRUE) res$proba <- "two-tailed probability" else res$proba <- "one-tailed probability"
	res$units.text <- UnitTxt
	attr(res, "units") <- Unit
	
	# Do we need to plot the graph for the turnogram?
	if (plotit == TRUE) {
		Ilevel <- -log(level, base=2)
		if (xlog == TRUE) xlogstr <- "x" else xlogstr <- ""
		if (two.tailed == TRUE) imin <- -1.1*Ilevel else imin <- 0
		subtext <- paste(fun, "/", res$proba)
		if (complete == FALSE) {
			yrange.dat <- c(res$info, imin, 1.1*Ilevel)
			yrange <- c(min(yrange.dat), max(yrange.dat))
			plot(res$interval, res$info, type="l", log=xlogstr, ylim=yrange, xlab=paste("interval", UnitTxt), ylab="I (bits)", main=paste("Simple turnogram for:", data), sub=subtext)
		} else {
			yrange <- c(min(c(res$info.min, imin)), max(c(res$info.max, 1.1*Ilevel)))
			plot(res$interval, res$info, type="l", log=xlogstr, ylim=yrange, xlab=paste("interval (", UnitTxt, ")", sep=""), ylab="I (bits)", main=paste("Complete turnogram for:", data), sub=subtext)
			lines(res$interval, res$info.min)
			lines(res$interval, res$info.max)
		}
		if (lhorz == TRUE) {
			if (two.tailed == TRUE) {
				abline(h=0)
				abline(h=-Ilevel, lty=2, col=2)
			}
			abline(h=Ilevel, lty=2, col=2)
		}
		if (lvert == TRUE) abline(v=Level, lty=2, col=2)
	}
	class(res) <- "turnogram"
	res 	# Return results
}
