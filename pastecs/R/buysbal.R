"buysbal" <-
function(x, y=NULL, frequency=12, units="years", datemin=NULL, dateformat="m/d/Y", count=FALSE) {
	# Check the frequency argument
	if (!is.numeric(frequency))
		stop("frequency must be a numeric value")
	frequency <- round(frequency)
	if (frequency < 2)
		stop("frequency must be a positive integer larger or equal to 2")
	# If y is not null, we should have two vectors: time in x, and data in y
	if (!is.null(y)) {
		if (!is.vector(x) || !is.vector(y) || length(x) != length(y))
			stop("x and y must be vectors of the same length")
		# Make sure x and y are vectors
		x <- as.vector(x)
		y <- as.vector(y)
	} else {		# We must have a time series in x
		if (!is.tseries(x))
			stop("x must be a regular time series if y is not provided")
		y <- as.vector(x)
		x <- as.vector(time(x))
	}
	# Verify units argument (currently, only "years" and "days" are accepted)
	UNITS <- c("years", "days")
	units <- pmatch(units, UNITS)
	if (is.na(units)) 
		stop("invalid units value")
	if (units == -1) 
		stop("ambiguous units value")
	x <- switch(units,
		"years"=x,									# No transformation required
		"days"=daystoyears(x, datemin, dateformat))	# Transform in years
	# Create factors for "years" and "cycle"
	years <- floor(x)
	cycle <- floor((x - years+0.00001)*frequency) + 1
	# Rem: we had to add a small value to cope with rounding errors with time series!
	# Check if all cycle values happen
	nocycle <- setdiff(1:frequency, cycle)
	# If at least one value does not happen, we must create dummy entries
	if (length(nocycle) > 0) {
		n <- length(nocycle)	
		cycle <- c(cycle, nocycle)
		years <- c(years, years[1:n])
		y <- c(y, rep(NA, n))
	}
	if (count == TRUE) { # Indicate the number of observations in each case
		y[!is.na(y)] <- 1
		tab <- tapply(y, list(years, cycle), FUN="sum", simplify=TRUE, na.rm=TRUE)
		tab[is.nan(tab)] <- 0
		tab[is.na(tab)] <- 0
	} else {			 # Create the Buys-Ballot table
		tab <- tapply(y, list(years, cycle), FUN="mean", simplify=TRUE, na.rm=TRUE)
		# We get mixed combinations of NA and NaN for missing entries => change NaN into NA
		tab[is.nan(tab)] <- NA
	}
	tab
}
