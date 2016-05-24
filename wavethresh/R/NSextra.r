"makewpstRO" <- 
function(timeseries, response, filter.number = 10., family = 
	"DaubExPhase", trans = logabs, percentage = 10.)
{
	#
	#
	# Using the data in time series (which should be a length a power of two)
	# and the response information. Create an object
	# of class wpstRO (stationary wavelet packet regression Object).
	#
	# Given this wpstRO and another timeseries a function exists to predict
	# the group membership of each timeseries element
	#
	#
	# First build stationary wavelet packet object
	#
	#
	# Now convert this to a matrix 
	#
	twpst <- #
	wpst(timeseries, filter.number = filter.number, family = family
		)
	#
	# Now extract the ``best'' 1D variables.
	#
	tw2m <- #
	wpst2m(wpstobj = twpst, trans = trans)
	tbm <- #
	bestm(tw2m, y = response, percentage = percentage)
	#
	# Now build data frame from these variables
	#
	#
	print.w2m(tbm)
	nc <- ncol(tbm$m)
	nr <- nrow(tbm$m)
	tdf <- data.frame(response, tbm$m)
	dimnames(tdf) <- list(as.character(1.:nr), c("response", paste(
		"X", 1.:nc, sep = "")))
	l <- list(df = tdf, ixvec = tbm$ixvec, level = tbm$level, pktix
		 = tbm$pktix, nlevels = tbm$nlevels, cv = tbm$cv, 
		filter = twpst$filter, trans = trans)
	oldClass(l) <- "wpstRO"
	l
}


"wpstREGR" <- 
function(newTS, wpstRO)
{
	#
	# Extract the "best packets"
	#
	newwpst <- #
	wpst(newTS, filter.number = wpstRO$filter$filter.number, family
		 = wpstRO$filter$family)
	goodlevel <- wpstRO$level
	goodpkt <- wpstRO$pkt
	npkts <- length(goodpkt)
	ndata <- length(newTS)
	m <- matrix(0., nrow = ndata, ncol = npkts)
	J <- nlevelsWT(newwpst)
	grot <- compgrot(J, filter.number=wpstRO$filter$filter.number,
		family=wpstRO$filter$family)
	for(i in 1.:npkts) {
		j <- goodlevel[i]
		m[, i] <- guyrot(accessD(newwpst, level = j, index = 
			goodpkt[i]), grot[J - j])/(sqrt(2.)^(J - j))
		m[, i] <- wpstRO$trans(m[, i])
	}
	dimnames(m) <- list(NULL, paste("X", 1.:npkts, sep = ""))
	l <- data.frame(m)
	l
}

"wpst2m" <- 
function(wpstobj, trans = identity)
{
	#
	# Function that converts a wpstobj into a matrix
	#
	# Input:
	#
	#		wpstobj:	the wpstobj to convert
	#
	#		trans:		the transform to apply to the
	#				wpst coefficients as they come out
	#
	#				an interesting alternative is
	#				trans = log( . )^2
	#				(you'll have to write this function)
	#
	#
	# Returns: An object of class w2m
	#
	#	This is a list with the following components:
	#
	#	m	- a matrix of order ndata x nbasis
	#
	#		where ndata is the number of data points for
	#		the time series that constituted wpstobj
	#
	#		and nbasis is the number of bases in the wpstobj 
	#
	#		Each column corresponds to a basis function
	#
	#		The row ordering is the same as the time series
	#		that constituted wpstobj
	#
	#	pktix	-	a vector of length nbasis which
	#			describes the packet index of the
	#			basis function in wpstm
	#
	#	level	-	as pktix but for the level
	#
	#	nlevels	The number of levels
	#
	J <- nlev <- nlevelsWT(wpstobj)
	grot <- compgrot(J, filter.number = wpstobj$filter$filter.number,
		family = wpstobj$filter$family)
	nbasis <- 2. * (2.^nlev - 1.)
	ndata <- 2.^nlev
	m <- matrix(0., nrow = ndata, ncol = nbasis)
	level <- rep(0., nbasis)
	pktix <- rep(0., nbasis)
	cnt <- 1.
	cat("Level: ")
	for(j in 0.:(nlev - 1.)) {
		cat(j, " ")
		lcnt <- 0.
		npkts <- 2.^(nlev - j)
		prcnt <- as.integer(npkts/10.)
		for(i in 0.:(npkts - 1.)) {
			pkcoef <- guyrot(accessD(wpstobj, level = j,
				index = i), grot[J - j])/(sqrt(2.)^
				(J - j))
			m[, cnt] <- trans(pkcoef)
			level[cnt] <- j
			pktix[cnt] <- i
			lcnt <- lcnt + 1.
			cnt <- cnt + 1.
			if (prcnt > 0)	{
				if(lcnt %% prcnt == 0.) {
					lcnt <- 0.
					cat(".")
					}
				}
		}
		cat("\n")
	}
	cat("\n")
	l <- list(m = m, level = level, pktix = pktix, nlevels = J)
	oldClass(l) <- "w2m"
	l
}

"compgrot" <- 
function(J, filter.number, family)
{
	if(filter.number == 1. && family == "DaubExPhase") {
		grot <- (2.^(0.:(J - 1.)) - 1.)
	}
	else {
		grot <- (1.:J)^2.
		grot[1.] <- 2.
		grot <- cumsum(grot)
	}
	grot
}

"logabs" <- 
function(x)
logb(x^2.)

"bestm" <- 
function(w2mobj, y, percentage = 50.)
{
	#
	# Compute desired number of bases
	#
	ndata <- #
	nrow(w2mobj$m)
	#
	# Actual number of bases
	#
	dbasis <- #
	as.integer((percentage * ndata)/100.)
	nbasis <- ncol(w2mobj$m)
	cv <- rep(0., nbasis)
	for(i in 1.:nbasis) {
		cv[i] <- cor(w2mobj$m[, i], y)
	}
	cv[is.na(cv)] <- 0.
	sv <- rev(sort.list(abs(cv)))[1.:dbasis]
	ixvec <- 1.:nbasis
	l <- list(m = w2mobj$m[, sv], ixvec = ixvec[sv], pktix = w2mobj$
		pktix[sv], level = w2mobj$level[sv], nlevels = w2mobj$
		nlevels, cv = cv[sv])
	oldClass(l) <- "w2m"
	l
}

"print.w2m" <- 
function(x, maxbasis = 10., ...)
{
	w2mobj <- x
	cat("Contains SWP coefficients\n")
	cat("Original time series length: ", nrow(w2mobj$m), "\n")
	cat("Number of bases: ", ncol(w2mobj$m), "\n")
	lbasis <- min(maxbasis, ncol(w2mobj$m))
	if(is.null(w2mobj$ixvec)) {
		cat("Raw object\n")
		mtmp <- cbind(w2mobj$level[1.:lbasis], w2mobj$pktix[
			1.:lbasis])
		dimnames(mtmp) <- list(NULL, c("Level", "Pkt Index"))
	}
	else {
		cat("Some basis selection performed\n")
		mtmp <- cbind(w2mobj$level[1.:lbasis], w2mobj$pktix[
			1.:lbasis], w2mobj$ixvec[1.:lbasis], w2mobj$
			cv[1.:lbasis])
		dimnames(mtmp) <- list(NULL, c("Level", "Pkt Index",
			"Orig Index", "Score"))
	}
	print(mtmp)
	if(ncol(w2mobj$m) > maxbasis)
		cat("etc. etc.\n")
	invisible()
}

"print.wpstRO" <- 
function(x, maxbasis = 10., ...)
{
	wpstRO <- x
	cat("Stationary wavelet packet regression object\n")
	cat("Composite object containing components:")
	cat("Original time series length: ", nrow(wpstRO$df), "\n")
	cat("Number of bases: ", ncol(wpstRO$df) - 1., "\n")
	lbasis <- min(maxbasis, ncol(wpstRO$df) - 1.)
	if(is.null(wpstRO$ixvec)) {
		cat("Raw object\n")
		mtmp <- cbind(wpstRO$level[1.:lbasis], wpstRO$pktix[
			1.:lbasis])
		dimnames(mtmp) <- list(NULL, c("Level", "Pkt Index"))
	}
	else {
		cat("Some basis selection performed\n")
		mtmp <- cbind(wpstRO$level[1.:lbasis], wpstRO$pktix[
			1.:lbasis], wpstRO$ixvec[1.:lbasis], wpstRO$
			cv[1.:lbasis])
		dimnames(mtmp) <- list(NULL, c("Level", "Pkt Index",
			"Orig Index", "Score"))
	}
	print(mtmp)
	if(ncol(wpstRO$df) > maxbasis)
		cat("etc. etc.\n")
	invisible()
}

