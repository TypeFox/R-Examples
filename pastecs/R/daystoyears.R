"daystoyears" <-
function (x, datemin=NULL, dateformat="m/d/Y") {
	x <- x
	datemin <- datemin
	defyearorig <- 1970
	# In R, we use POSIXt
	if (length(datemin) > 0 && !any(class(datemin) == "POSIXt")) {	# We must convert
		# To be compatible with chron() and with Splus,  we accept format as "d/m/y"
		# which is converted into %d/%m/%y. Warning! Still must make difference between
		# "y" that corresponds to year in "89" and "Y" for full-spelled years in "1989"!!!
		# make necessary conversions to accept old formats
		dateformat <- sub("month", "%B", dateformat)	# Full spelled month
		dateformat <- sub("mon", "%b", dateformat)		# Abbreviated month
		dateformat <- sub("m", "%m", dateformat)		# month - numeric
		dateformat <- sub("d", "%d", dateformat)	    # day
		dateformat <- sub("y", "%y", dateformat)		# year (two digits)
		dateformat <- sub("Y", "%Y", dateformat)		# Year (four digits)
		datemin <- strptime(as.character(datemin), format=dateformat)
	}
	# Julian is adapted from julian.default in lib chron 2.2-19 (this way we don't require chron!)
	"Julian" <- function(x, d, y) {
		if(is.null(origin. <- getOption("chron.origin")))
				origin. <- c(month = 1, day = 1, year = 1970)	# Default origin in R
		m <- c(origin.[1], x)               # prepend month of new origin
		d <- c(origin.[2], d)               # prepend day of new origin
		y <- c(origin.[3], y)               # prepend year of new origin
		# code from julian date in the S book (p.269)
		y <- y + ifelse(m > 2, 0, -1)
		m <- m + ifelse(m > 2, -3, 9)
		c <- y %/% 100
		ya <- y - 100 * c
		out <- ((146097 * c) %/% 4 + (1461 * ya) %/% 4 + (153 * m + 2) %/% 5 + d + 1721119)
		## now subtract the new origin from all dates
		if(all(origin. == 0))
			out <- out[-1]
		else
			out <- out[-1] - out[1]
		# orig according to S algorithm
		out
	}

	if (length(datemin) > 0) {
		dateminval <-  Julian(datemin$mon+1, datemin$mday, datemin$year+1900)
		# now we shift the whole x series so as the minimal day matches dateminval
		x <- x - trunc(min(x, na.rm=TRUE)) + dateminval
	}
	# We have days as units. We want years with a "linear scale", i.e.: 1 year = 365.25 days, 1 month = 1/12 years
	# We want also the integer value reflect exactly the current year, i.e.: 1997.xxx for dates in the year 1997
	if(is.null(yearorig <- options("chron.origin")$year))
		yearorig <- defyearorig	
	x <- x/365.25 + yearorig
	x
}
