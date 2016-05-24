"regul" <-
function (x, y = NULL, xmin=min(x), n=length(x), units="days", frequency=NULL, deltat=1/frequency, datemin=NULL, dateformat="m/d/Y", tol=NULL, tol.type="both", methods="linear", rule=1, f=0, periodic=FALSE, window=(max(x)-min(x))/(n-1), split=100, specs=NULL) {
	call <- match.call()
	# Do we have specs?
	if (!is.null(specs)) {
		# Verify it is an objct of the right class 'specs.regul'
		specs.class <- class(specs)
		if (is.null(specs.class) || specs.class != "specs.regul")
			stop("specs must be a 'specs.regul' object.")
		# For each argument we look if it was not explicitly given (in this case, we take value in specs)
		arg.names <- names(call)
		if (pmatch("xmin", arg.names, 0) == 0)
			xmin <- specs$xmin
		if (pmatch("n", arg.names, 0) == 0)
			n <- specs$n
		if (pmatch("frequency", arg.names, 0) == 0)
			frequency <- specs$frequency
		if (pmatch("deltat", arg.names, 0) == 0)
			deltat <- specs$deltat
		if (pmatch("units", arg.names, 0) == 0)
			units <- specs$units
		if (pmatch("datemin", arg.names, 0) == 0)
			datemin <- specs$datemin
		if (pmatch("dateformat", arg.names, 0) == 0)
			dateformat <- specs$dateformat
		if (pmatch("tol", arg.names, 0) == 0)
			tol <- specs$tol
		if (pmatch("tol.type", arg.names, 0) == 0)
			tol.type <- specs$tol.type
		if (pmatch("methods", arg.names, 0) == 0)
			methods <- specs$methods
		if (pmatch("rule", arg.names, 0) == 0)
			rule <- specs$rule
		if (pmatch("f", arg.names, 0) == 0)
			f <- specs$f
		if (pmatch("periodic", arg.names, 0) == 0)
			periodic <- specs$periodic
		if (pmatch("window", arg.names, 0) == 0)
			window <- specs$window
		if (pmatch("split", arg.names, 0) == 0)
			split <- specs$split	
	}
	# Evaluate arguments now
	x <- x
	y <- y
	xmin <- xmin
	n <- n
	units <- units
	frequency <- frequency
	deltat <- deltat
	# tol <- tol		# Must be evaluated only after recalculation of deltat!
	tol.type <- tol.type
	datemin <- datemin
	dateformat <- dateformat
	methods <- methods
	rule <- rule
	f <- f
	periodic <- periodic
	window <- window
	split <- split
	# Verify arguments that are not verified further down
	# x, y, xmin, tol.type, methods, periodic, window, split are verified elsewhere
	# we don't care for units and rule
	# Set frequency and deltat correctly
	if (is.null(frequency)) {
		if (is.null(deltat))
			deltat <- 1					# value by default when none is provided = 1
			frequency <- 1/deltat
		} else {							# frequency provided. Make sure deltat is its inverse
			deltat <- 1/frequency
	}
	# deltat must be a strictly positive value
	if (deltat <= 0)
		stop("deltat must be a positive number")
	# Make sure tol is defined
	tol <- tol
	if (is.null(tol)) tol <- 0
	# n must be a positive integer
	if (n <= 0)
		stop("n must be a positive number")
	if (n != round(n))
		stop("n must be an integer number")
	# create specs
	specs <- list(xmin=xmin, n=n, frequency=frequency, deltat=deltat, units=units, datemin=datemin, dateformat=dateformat, tol=tol, tol.type=tol.type, methods=methods, rule=rule, f=f, periodic=periodic, window=window, split=split)
	# Make sure data are sorted in increasing order with time
	srt <- sort.list(x)
	x <- x[srt]
	if (is.null(ncol(y))) y <- y[srt] else y <- y[srt,]
	# If datemin is provided, we shift x accordingly
	if (length(datemin)>0) {
		dateminval <- yearstodays(daystoyears(1, datemin, dateformat))
		xminval <- trunc(min(x, na.rm=TRUE))
		x <- x - xminval + dateminval
		xmin <- xmin - xminval + dateminval
	}
	# If units is "daystoyears", there is a special treatment to transform the time-scale!
	if (units == "daystoyears") {
		# We have days as units. We want years with a "linear scale", i.e.: 1 year = 365.25 days, 1 month = 1/12 years
		# We want also the integer value reflect exactly the current year, i.e.: 1997.xxx for dates in the year 1997
		if(is.null(yearorig <- options("chron.origin")$year)) {
			yearorig <- 1970
		}
		x <- x/365.25 + yearorig
		xmin <- xmin/365.25 + yearorig
		cat(paste("A 'tol' of", tol, "in 'days' is", tol/365.25, "in 'years'\n"))
		tol <- tol/365.25
		window <- window/365.25
		units <- "years"					# and we reset units to "years"
	}
	# Verify the methods argument
	METHODS <- c("constant", "linear", "spline", "area")
	# repeat methods a certain times to make sure it match at least the number of columns to regulate
	nser <- ncol(y)
	if (is.null(nser)) nser <- 1		# when y is a vector
	if (length(methods) < nser)
		methods <- rep(methods, length.out=nser)	
	methindex <- NULL
	for (i in 1:length(methods)) {
		methindex[i] <- pmatch(methods[i], METHODS)
		if (is.na(methindex[i])) 
			stop(paste("invalid regulation method:", methods[i]))
		if (methindex[i] == -1) 
			stop(paste("ambiguous regulation method:", methods[i]))
	}
	# Regulate each column of the matrix in turn
	# If y is a vector, only one call is required
	if (is.vector(y)) {
		# Choose the method
		res <- switch(methindex[1],
		    	"constant"=regconst(x, y, xmin=xmin, n=n, deltat=deltat, rule=rule, f=f),
		    	"linear"=reglin(x, y, xmin=xmin, n=n, deltat=deltat, rule=rule),
				"spline"=regspline(x, y, xmin=xmin, n=n, deltat=deltat, rule=rule, periodic=periodic),
				"area"=regarea(x, y, xmin=xmin, n=n, deltat=deltat, rule=rule, window=window, interp=TRUE, split=split))
		xout <- res$x
		yout <- res$y
		yout <- as.data.frame(yout)
		names(yout) <- "Series"
	} else {
		# If x is not a vector, it is treated as a data frame
		# A result will be returned in a data frame with corresponding columns
		y <- as.data.frame(y)
		# We keep the same column headers
		NamesV <- names(y)
		yout <- NULL
		# Calculation is performed alternatively on each column
		for (i in 1:nser) {
			res <- switch(methindex[i],
		    		"constant"=regconst(x, y[,i], xmin=xmin, n=n, deltat=deltat, rule=rule, f=f),
		    		"linear"=reglin(x, y[,i], xmin=xmin, n=n, deltat=deltat, rule=rule),
					"spline"=regspline(x, y[,i], xmin=xmin, n=n, deltat=deltat, rule=rule, periodic=periodic),
					"area"=regarea(x, y[,i], xmin=xmin, n=n, deltat=deltat, rule=rule, window=window, interp=TRUE, split=split))
			# The next if condition to avoid error at the first deltat
			if (is.null(yout)==TRUE) yout <- data.frame(res$y) else yout <- cbind(yout, res$y)
		}
		# Change names of columns to match the original data frame
		names(yout) <- NamesV
		# Get xout
		xout <- res$x
	}
	# make sure methods are fully spelled
	methods <- NULL
	for (i in 1:nser) {
		methods[i] <- switch(methindex[i],
			    	"constant"="constant",
			    	"linear"="linear",
					"spline"="spline",
					"area"="area")
	}
	# make sure yini is a date frame
	yini <- as.data.frame(y)
	if (ncol(yini) == 1) names(yini) <- "Initial"
	# now we calculate match according to tol.type
	TOL.TYPES <- c("left", "both", "right", "none")
	tol.idx <- pmatch(tol.type, TOL.TYPES)
	if (is.na(tol.idx)) 
		stop("invalid tol.type value")
	if (tol.idx == -1) 
		stop("ambiguous tol.type value")
	# make sure tol.type is fully spelled
	tol.type <- switch(tol.idx,
			    "left"="left",
			    "both"="both",
				"right"="right",
				"none"="none")
	if (is.null(tol)) tol <- 0 else tol <- abs(tol)
	pos <- match.tol(xout, x, tol.type=tol.type, tol=tol)
	# We must recalculate tol here, to keep the right value
	if (tol > deltat) tol2 <- deltat else tol2 <- deltat/round(deltat/tol)
	# we replace corresponding values for each series in yout by indexed values of the same series in y
	repl <- yini[pos, ]
	yout <- as.matrix(yout)			# The next instruction does not work on data frames in R!
	yout[!is.na(repl)] <- repl[!is.na(repl)]
	yout <- as.data.frame(yout)
	# which regular date match observations? Put distance for a matching observation, Inf for an interpolated observation, -Inf for an extrapolated observation
	match.dist <- xout - x[pos]
	match.dist[is.na(match.dist)] <- Inf
	match.dist[(xout < min(x, na.rm=TRUE) | xout > max(x, na.rm=TRUE)) & match.dist == Inf] <- -Inf
	# create tspar
	# start is a numeric vector of 2 numbers indicating integer and fractional part of start.
	# its combination must be a multiple of deltat => adjust to the closest neighbour
	dateini <- xout[1]
	intini <- trunc(dateini)									# Integer part of start
	# decini must be lower than frequency
	if (intini != dateini & frequency > 1) {
		decini <- floor((dateini - intini) * frequency) + 1			# "Decimal" part of start must be a number times of frequency for 'ts' objects
	} else {
		decini <- NULL
	}
	tspar <- list("start"=c(intini, decini), "deltat"=deltat, "frequency"=frequency)
	# create the list containing the results
	res <- list(x=xout, y=yout, xini=x, yini=yini, tspar=tspar, units=units, match.dist=match.dist, specs=specs, call=call)
	class(res) <- "regul"										# and turn it into a 'regul' object
	res															# Return the result
}
