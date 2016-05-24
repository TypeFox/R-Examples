"match.tol" <-
function(x, table, nomatch=NA, tol.type="both", tol=0) {
	# replace interpolated/extrapolated value yout(i) by a real observation if there is one in the time scale
	# xout(i) +/- tol if tol.type = "center"
	# xout(i) + tol if tol.type = "right"
	# xout(i) - tol if tol.type = "left"
	# nothing if tol.type = "none"
	# and by an exact matching if tol = 0 and tol.type != "none"
	# If, after rounding, there are several candidates for matching, only the closest value is kept
	# The matrix algorithm used here work only for round fraction of deltat for tol!! (tol = deltat, tol = deltat/2, etc...)
	
	# Right match (in case of several matches, the closest is the first one, reported correctly by match)
	match.tol.right <- function(x, table, tol, nomatch=NA) {
		if (tol == 0) {
			posr <- match(x, table, nomatch)
		} else {			# Warning: tol must be a round fraction of deltat
			xornd <- round((x - min(x, na.rm=TRUE)) / tol)
			xrnd <- floor((table - min(x, na.rm=TRUE)) / tol)
			posr <- match(xornd, xrnd, nomatch)
		}
		posr
	}
	# Left match (in case of several matches, the closest is the last one => apply match on the reversed vector)
	match.tol.left <- function(x, table, tol, nomatch=NA) {
		if (tol == 0) {
			posl <- match(x, table, nomatch)
		} else {			# Warning: tol must be a round fraction of deltat
			xornd <- round((x - min(x, na.rm=TRUE)) / tol)
			xrnd <- ceiling((table - min(x, na.rm=TRUE)) / tol)
			xol <- length(xornd)
			xl <- length(xrnd)
			posl <- match(xornd[xol:1], xrnd[xl:1], nomatch)
			# posl is reversed and indices are changed to correspond to a usual match, but keeping last one in case of ex aequos
			posl <- xl + 1 - posl[xol:1]
		}
		posl
	}
	# match on both sides. We search for tol on left and on right
	match.tol.both <- function(x, table, tol, nomatch=NA) {
		if (tol == 0) {
			posr <- match(x, table, nomatch)
		} else {			# Warning: tol must be a round fraction of deltat
			# Seach right
			xornd <- round((x - min(x, na.rm=TRUE)) / tol)
			xrnd <- floor((table - min(x, na.rm=TRUE)) / tol)
			posr <- match(xornd, xrnd, nomatch)
			# Search left
			xrnd <- ceiling((table - min(x, na.rm=TRUE)) / tol)
			xol <- length(xornd)
			xl <- length(xrnd)
			posl <- match(xornd[xol:1], xrnd[xl:1], nomatch)
			# posl is reversed and indices are changed to correspond to a usual match, but keeping last one in case of ex aequos
			posl <- xl + 1 - posl[xol:1]
			# we keep only posr, but replace NA by values of posl, and also in case of both left, and right matches
			# we compare which one is closest. If left match is closest than right match, we replace also in posr
			repl <- (is.na(posr)) | ((x - table[posl]) < (table[posr] - x))		# In case of same distance, we keep right match
			repl[is.na(repl)] <- FALSE
			posr[repl] <- posl[repl]
		}
		posr
	}
				
	# match.tol starts here
	table <- sort(table)					# Make sure table is sorted in increasing order
	table <- table[!is.na(table)]			# Eliminate missing values
	n <- length(table)
	# does x be regularly spaced?
	spaces <- x[2:length(x)] - x[1:(length(x)-1)]
	if (max(spaces) - min(spaces) > max(spaces)/100)
		stop("x must be a regularly spaced vector for match.tol!")
	# We verify also that tol is a round fraction of the space in x
	space <- mean(spaces)
	if (is.null(tol) || tol == 0) {
		tol <- 0
		tol2 <- 0
	} else {
		tol2 <- abs(tol)
		if (tol2 > space) tol2 <- space else {
			tol2 <- space/round(space/tol2)
		}
	}
	if (tol2 != tol)
		cat(paste("'tol' was adjusted to", tol2, "\n\n"))
	# match according to tol.type
	TOL.TYPES <- c("left", "both", "right", "none")
	tol.idx <- pmatch(tol.type, TOL.TYPES)
	if (is.na(tol.idx)) 
		stop("invalid tol.type value")
	if (tol.idx == -1) 
		stop("ambiguous tol.type value")
	# calculate the matching vector
	pos <- switch(tol.idx,
		  	"left"=match.tol.left(x, table, tol2),
		  	"both"=match.tol.both(x, table, tol2),
			"right"=match.tol.right(x, table, tol2),
			"none"=match(x, table, nomatch=NA))
	pos
}
