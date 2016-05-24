"extract.escouf" <-
function(e, n, level=e$level, ...) {
	if (missing(n)) n <- NULL
	if (!is.null(n)) {		# We want to extract n variables
		# We calculate the corresponding level
		n <- abs(round(n))			# Make sure n is a positive integer!
		if (n>length(e$RV)-1) {		# We extract all variables
			level <- 1
		} else {					# We calculate the level
			level <-(e$RV[n]+e$RV[n+1])/2
		}
	}
	if (is.null(level)) {		# look if object$level exist
		if (is.null(e$level)) {
			stop("You must provide a level value for extraction!")
		} else {
			level <- e$level
		}
	}
	# Check the validity of level
	if (level>1 || level<0) stop("level must be a value between 0 and 1!")
	# level must not be lower than e$calc.level, otherwise we don't have enough information!
	if (level>e$calc.level) stop("level is higher that the one used in calculation, unable to fully extract Escoufier's matrix at this level!")
	Res <- eval(parse(text=e$data))[e$vr[e$RV<level]]
	Res
}
