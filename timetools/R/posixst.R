POSIXst <- function (x, unit, of=NULL, tz='UTC', ...) UseMethod('POSIXst')

POSIXst.default <- function (x, unit, of=NULL, tz='UTC', ...)
{
	if (missing (x))
		return(POSIXst (as.POSIXlt(character(), tz),
				unit, of, tz, ...))
	stop(sprintf("'POSIXst' method not implemented for %s object", class(x)))
}

POSIXst.integer <- function (x, unit, of=NULL, tz='UTC', ...)
{
	unit <- POSIXt.units(unit)
	if( unit == POSIXt.units('year') & is.null(of) )
		of <- POSIXt.units('AD') else
	if( unit == POSIXt.units('month') & is.null(of) )
		of <- POSIXt.units('year') else
		of <- POSIXt.units(of)
	new('POSIXst', subtime=x, unit=unit, of=of, timezone=tz)
}

POSIXst.numeric <- function (x, unit, of=NULL, tz='UTC', ...)
{
	x <- as.character(x)
	if( any(grepl('\\.', x)) ) stop("'x', must be a whole number)")
	x <- as.integer(x)
	POSIXst(x, unit, of, tz, ...)
}

POSIXst.POSIXct <- function (x, unit, of=NULL, tz=attributes(x)$tzone, ...)
{
	x <- as.POSIXlt(x)
	POSIXst (x, unit, of, tz[1], ...)
}

POSIXst.POSIXlt <- function (x, unit, of=NULL, tz=attributes(x)$tzone, ...)
{
	unit <- POSIXt.units(unit)
	if( unit == POSIXt.units('year') & is.null(of) )
		of <- POSIXt.units('AD') else
	if( unit == POSIXt.units('month') & is.null(of) )
		of <- POSIXt.units('year') else
		of <- POSIXt.units(of)

	x <- as.POSIXlt(as.POSIXct(x), tz[1])

	if (unit == POSIXt.units('year')) {
		res <- x$year+1900
	} else if (unit==POSIXt.units('month')) {
		res <- x$mon
	} else if (unit == POSIXt.units('day')) {
		res <- switch (as.character(of),
			year = x$yday,
			month = x$mday,
			week = x$wday,
			"'of' should be one of (year, month, week)")
	} else if (unit == POSIXt.units('hour')) {
		res <- switch (as.character(of),
			year	=  x$yday	* 24 + x$hour,
			month	= (x$mday-1)	* 24 + x$hour,
			week	=  x$wday	* 24 + x$hour,
			day	= 		       x$hour,
			"'of' should be one of (year, month, week, day)")
	} else if (unit == POSIXt.units('minute')) {
		res <- switch (as.character(of),
			year	=( x$yday	* 24 + x$hour) * 60 + x$min,
			month	=((x$mday-1)	* 24 + x$hour) * 60 + x$min,
			week	=( x$wday	* 24 + x$hour) * 60 + x$min,
			day	= (		       x$hour) * 60 + x$min,
			hour	=				      x$min,
			"'of' should be one of (year, month, week, day, hour)")
	} else if (unit == POSIXt.units('second')) {
		res <- switch (as.character(of),
			year	=( (x$yday	* 24 + x$hour) * 60 + x$min) * 60 + x$sec,
			month	=(((x$mday-1)	* 24 + x$hour) * 60 + x$min) * 60 + x$sec,
			week	=( (x$wday	* 24 + x$hour) * 60 + x$min) * 60 + x$sec,
			day	=((		       x$hour) * 60 + x$min) * 60 + x$sec,
			hour	=(				      x$min) * 60 + x$sec,
			minute	= 						    x$sec,
		"'of' should be one of (year, month, week, day, hour, minute)")
	}
	return( POSIXst(res, unit, of, tz) )
}

POSIXst.TimeInstantDataFrame <-
	function(x, unit, of=NULL, tz=timezone(x), ...) {
	POSIXst(force(when(x)), unit, of, tz, ...)
}

POSIXst.TimeIntervalDataFrame <-
	function(x, unit, of=NULL, tz=timezone(x), ..., cursor=NULL) {
	if( !is.null(cursor) )
		return( POSIXst(as.TimeInstantDataFrame(x, cursor),
				unit, of, tz) )

	unit <- POSIXt.units(unit)
	u <- as.character(unit)
	u <- switch(u,
		second='sec', minute='min', u)
	st <- mapply(seq.POSIXt, SIMPLIFY=FALSE,
		     start(x), end(x),
		     MoreArgs=list(by=u))

	st <- lapply(st, POSIXst, unit, of, tz)

	return( st )
}

# intuitive wrappers

year	<- function(x, ...) UseMethod('year')
month	<- function(x, ...) UseMethod('month')
day	<- function(x, of, ...) UseMethod('day')
hour	<- function(x, of, ...) UseMethod('hour')
minute	<- function(x, of, ...) UseMethod('minute')
second	<- function(x, of, ...) UseMethod('second')

setGeneric ('year',	function (x, ...) standardGeneric ('year') )
setGeneric ('month',	function (x, ...) standardGeneric ('month') )
setGeneric ('day',	function (x, of, ...) standardGeneric ('day') )
setGeneric ('hour',	function (x, of, ...) standardGeneric ('hour') )
setGeneric ('minute',	function (x, of, ...) standardGeneric ('minute') )
setGeneric ('second',	function (x, of, ...) standardGeneric ('second') )

setMethod ('year', 'ANY',
	   function(x, ...) POSIXst(x, unit='year', ...) )
setMethod ('month', 'ANY',
	   function(x, ...) POSIXst(x, unit='month', ...) )
setMethod ('day', 'ANY',
	   function(x, of, ...) POSIXst(x, unit='day', of=of, ...))
setMethod ('hour', 'ANY',
	   function(x, of, ...) POSIXst(x, unit='hour', of=of, ...))
setMethod ('minute', 'ANY',
	   function(x, of, ...) POSIXst(x, unit='minute', of=of, ...))
setMethod ('second', 'ANY',
	   function(x, of, ...) POSIXst(x, unit='second', of=of, ...))

# TODO : define 'timezone' and 'timezone<-' for POSIXct.
