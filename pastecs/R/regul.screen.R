"regul.screen" <-
function(x, weight=NULL, xmin=min(x), frequency=NULL, deltat=(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))/(length(x)-1), tol=deltat/5, tol.type="both") {

	regul.screen.calc <- function(x, weight, xmin, deltat, tol, tol.type) {
		xmin <- xmin
		deltat <- deltat
		# We verify that tol is a round fraction of deltat
		if (is.null(tol) || tol == 0) tol2 <- 0 else {
			tol2 <- abs(tol)
			if (tol2 > deltat) tol2 <- deltat else {
				tol2 <- deltat/round(deltat/tol2)
			}
		}
		# We calculate highest n, so as there is no extrapolation
		if (max(x) < xmin) {		# There is no solution!
			nx <- 0
			stop("xmin is higher than max value in x!")
		} else {					# We calculate nx
			nx <- floor((max(x) - xmin)/deltat) + 1
		}
		# create xout vector
		xout <- 0:(nx-1) * deltat + xmin
		# calculate the matching vector
		pos <- match.tol(xout, x, tol.type=tol.type, tol=tol2)
		# which regular date match observations? Put distance for a matching observation, Inf for an interpolated observation, -Inf for an extrapolated observation
		match.dist <- xout - x[pos]
		match.dist[is.na(match.dist)] <- Inf
		match.dist[(xout < min(x) | xout > max(x)) & match.dist == Inf] <- -Inf
		# Weighted number of matches
		weight.vec <- weight[pos]
		match.vec <- as.numeric(is.finite(match.dist))
		match <- sum(match.vec * weight.vec, na.rm=TRUE)
		exact.match.vec <- as.numeric(match.dist == 0)
		exact.match <- sum(exact.match.vec * weight.vec,na.rm=TRUE)
		res <- list(tol=tol2, nx=nx, match=match, exact.match=exact.match)
		res
	}
	
	# regul.screen starts here
	if (!is.numeric(x))
		stop ("x must be a numerical vector")
	if (length(x) < 2)
		stop ("x must contain at least two values")
	keep <- !is.na(x)			# Eliminate missing values
	x <- x[keep]
	if (length(x) < 2)
		stop ("x must contain at least two non missing values")
	if (is.null(weight)) {										# create a uniform weight vector of 1's
		weight <- rep(1, length.out=length(x))
	} else {													# eliminate values corresponding to missing data in x
		weight <- weight[keep]
	}
	weight[is.na(weight)] <- 0									# replace missing values in weight by 0 (corresp. value in x will be ignored)
	weight[weight < 0] <- 0										# do the same for negative values
	if (length(weight) != length(x))
		stop("x and weight must be vectors of equal length")
	# make sure x is sorted in increasing order
	srt <- sort.list(x)
	x <- x[srt]
	weight <- weight[srt]
	if (is.null(deltat)) {
		if (is.null(frequency)) {
			stop("You must define at least one of frequency or deltat")
		} else {
			deltat <- 1/frequency
		}
	}
	if (is.null(tol) | sum(is.na(tol)) > 0) tol <- deltat/5		# Default value
	# We need a tol value for each deltat
	tol <- rep(tol, length.out=length(deltat))
	# We loop on the various values of deltat, and then on the different values of xmin
	# and report three tables: 1) nx, 2) total.matches, 3) total.exact.matches
	nx <- length(xmin)
	nd <- length(deltat)
	rnames <- paste("x=", xmin, sep="")
	cnames <- paste("d=", deltat, sep="")
	dnames <- list(rnames, cnames)
	ttol <- NULL
	tnx <- matrix(nrow=nx, ncol=nd, dimnames=dnames)
	ttotal <- matrix(nrow=nx, ncol=nd, dimnames=dnames)
	texact <- matrix(nrow=nx, ncol=nd, dimnames=dnames)
	for (i in 1:nx) {
		for (j in 1:nd) {
			res <- regul.screen.calc(x, weight, xmin=xmin[i], deltat=deltat[j], tol=tol[j], tol.type=tol.type)
			ttol[j] <- res$tol
			tnx[i,j] <- res$nx
			ttotal[i,j] <- res$match
			texact[i,j] <- res$exact.match
		}
	}
	names(ttol) <- cnames
	res <- list(tol=ttol, n=tnx, nbr.match=ttotal, nbr.exact.match=texact)
	res
}
