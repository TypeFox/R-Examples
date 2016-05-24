"decaverage" <-
function(x, type="additive", order=1, times=1, sides=2, ends="fill", weights=NULL) {
	call <- match.call()
	x <- as.ts(x)
	if (is.matrix(x) && ncol(x) != 1)
	    stop("only univariate series are allowed")
	if (!is.numeric(times) || times <= 0)
		stop("times must be a positive number")
	if (!is.numeric(sides) || (sides != 1 & sides != 2))
		stop("specify only 1 or 2 for sides")
	# if weights exists, we use it in priority. Otherwise, we build the weights specification
	if (is.null(weights)) {
		if (is.character(order)) {
			if (is.na(pmatch(order, "periodic")))
				stop("order must be a positive number or \"periodic\"")
			freq <- frequency (x)
			# frequency must be at least 2
			if (freq < 2)
				stop("for a periodic smoothing, frequency of the series must be equal or higher than 2")
			order <- freq %/% 2
			weights <- rep(1, 2*order+1)
			if (freq == 2*order) {		# freq is even => weights is odd, and we must give 0.5 for weight at each end
				weights[1] <- 0.5
				weights[length(weights)] <- 0.5
			}
		} else {		# Order is numeric, so it must be a positive number
			if (!is.numeric(order) || order <= 0)
				stop("order must be a positive number or \"periodic\"")
			weights <- rep(1,2*order+1)
		}
	} else {			# weights is defined
		if (length(weights) < 2)
			stop("weights must contain at least 2 elements")
		order <- length(weights) %/% 2
	}
	# Check the type argument
	TYPES <- c("additive", "multiplicative")
		typeindex <- pmatch(type, TYPES)
		if (is.na(typeindex)) 
			stop("invalid type value")
		if (typeindex == -1) 
			stop("ambiguous type value")
		# make sure type is fully spelled
		type <- switch(typeindex,
				"additive"="additive",
				"multiplicative"="multiplicative")
	# Check the ends argument and treat the series accordingly (add calculated arguments at the beginning and at the end)
	ENDS <- c("NAs", "fill", "circular", "periodic")
	endsindex <- pmatch(ends, ENDS)
	if (is.na(endsindex)) 
		stop("invalid ends value")
	if (endsindex == -1) 
		stop("ambiguous ends value")
	# make sure ends is fully spelled
	ends <- switch(endsindex,
				"NAs"="NAs",
				"fill"="fill",
				"circular"="circular",
				"periodic"="periodic")
	# create our own specs component
	specs <- list(method="average", type=type, order=order, times=times, sides=sides, ends=ends, weights=weights)
	# we recuperate units from x
	units <- attr(x, "units")
	# We define functions that pads elements at end of the vector
	padmean <- function(x, order, sides) {
		n <- length(x)
		if (sides == 2)	{		# pads at each end
			if (order > n) order <- n
			paddedx <- NULL
			paddedx[(1:n)+order] <- x
			paddedx[1:order] <- mean(x[1:order], na.rm=TRUE)
			paddedx[(1:order)+n+order] <- mean(x[(n - order + 1):n], na.rm=TRUE)
			# Rem: we dont change tspar, because we will eliminate these values latter!
			cut <- c(order+1, order+n)
		} else {				# pads only at left
			if (2*order > n) order <- n/2
			paddedx <- NULL
			paddedx[(1:n)+2*order] <- x
			paddedx[1:(2*order)] <- mean(x[1:(2*order)], na.rm=TRUE)
			# Rem: we dont change tspar, because we will eliminate these values latter!
			cut <- c(2*order+1, 2*order+n)
		}
		res <- list(x=paddedx, circular=FALSE, cut=cut)
		res
	}
	# In the next function, we take the equivalent sequence of first year to pad before first year
	# and the equivalent function of last year to pad after last year
	padper <- function(x, sides) {
		n <- length(x)
		f <- frequency(x)
		# We must have at least two complete cycles here!
		if (n < 2*f) {
			warning("you need at least two complete cycles to use ends = \"periodic\"")
			# We don't change the series
			res <- list(x=x, circular=FALSE, cut=c(1,n))
		} else {		# Requirements to calculate padper are met
			if (sides == 2)	{		# pads at each end
				pos <- f %/% 2
				pos0 <- f - pos + 1
				paddedx <- NULL
				paddedx[(1:n)+pos] <- x
				paddedx[1:pos] <- x[pos0:f]
				paddedx[(1:pos)+n+pos] <- x[(1:pos)+n-f]
				# Rem: we dont change tspar, because we will eliminate these values later!
				cut <- c(pos+1, pos+n)
			} else {				# pads only at left
				pos <- (f %/% 2) * 2
				pos0 <- f - pos + 1
				paddedx <- NULL
				paddedx[(1:n)+pos] <- x
				paddedx[1:pos] <- x[pos0:f]
				# Rem: we dont change tspar, because we will eliminate these values later!
				cut <- c(pos+1, pos+n)
			}
			res <- list(x=paddedx, circular=FALSE, cut=cut)
		}
		res
	}
	n <- length(x)
	filtered <- x						# We don't change the initial series, but a copy of it
	filt <- weights/sum(weights)		# Scale down weights
	for (i in 1:times) {
		# Pad elements at ends of vector x
		padx <- switch(endsindex,
				"NAs"=list(x=filtered, circular=FALSE, cut=c(1,n)),		# We don't have to change the series
				"fill"=padmean(filtered, order, sides),					# We take the mean of ends elements and put them at extremes
				"circular"=list(x=filtered, circular=TRUE, cut=c(1,n)),	# We don't change the series, but change circular (done latter)
				"periodic"=padper(filtered, sides))						# We add elements from first and last period at the right place
		circular <- padx$circular
		cut <- padx$cut
		# perform filtering
		filtered <- filter(padx$x, filter=filt, method="convolution", sides=sides, circular=circular)
		# Now we have to cut the vector x according to cut (we don't use the function window for that since we didn't changed tspar!)
		filtered <- as.ts(as.vector(filtered)[cut[1]:cut[2]])
		tsp(filtered) <- tsp(x)
		filtered<- as.ts(filtered)
	}
	# Calculate residuals
	if (type == "additive") {
		residuals <- x - filtered
	} else {
		residuals <- x / filtered
	}
	series <- ts.union(filtered, residuals)
	# create our own 'tsd' structure
	res <- list(ts="series", series=series, weights=weights, units=units, specs=specs, call=call)
	class(res) <- "tsd"		# change the class of the object to 'tsd'
	res
}
