# definition de la classe
#------------------------
setClass ('POSIXctp',
	representation (duration='integer', unit='factor'),
	validity=function(object) {
		if (length (object@unit) != length (object@duration) ) {
			return ("Slots 'duration' and 'unit' must have same length.")
		} else if (!all (object@unit %in% POSIXt.units() ) ) {
			return (sprintf ('%s are not a valid units.', paste (as.character (object@unit), collapse=', ') ) )
		} else return (TRUE)
	})

# constructeur
#-------------
POSIXctp <- function (duration, unit)
{
	if( missing (duration) ) duration <- 1L
	if( missing(unit) ) {
		unit <- duration
		duration <- rep(1L, length(unit))
	}

	if (inherits(duration, 'numeric') | inherits(duration, 'logical')) {
		duration <- as.integer (duration)
	}

	if( length(duration) > 1 & length(unit) == 1) {
		unit <- rep( unit, length(duration) )
	} else if( length(duration) < length(unit) ) {
		warning("'duration' is shorter than 'unit'. 'duration' is recycled.")
		duration <- rep(duration, length.out=length(unit))
	} else if( length(duration) > length(unit) ) {
		warning("'unit' is shorter than 'duration'. 'unit' is recycled.")
		unit <- rep(unit, length.out=length(duration))
	}

	if (inherits (unit, 'character'))
		unit <- POSIXt.units(unit)

	return (new('POSIXctp', duration=duration, unit=unit))
}

# acces aux proprietes
#---------------------
setMethod ('unit', 'POSIXctp', function(x, ...) (x@unit))

setMethod (f='unit<-', signature='POSIXctp',
definition=function(object, value) {
	value <- POSIXt.units(value)
	conversion <- array (
		c(1, 12, rep(NA, 5),
		  NA, 1, rep(NA, 5), 
		  rep(NA, 2), 1, 7, 7*24, 7*24*60, 7*24*60*60,
		  rep(NA, 3), 1, 24, 24*60, 24*60*60,
		  rep(NA, 4), 1, 60, 60*60,
		  rep(NA, 5), 1, 60, rep(NA, 6), 1),
		dim=c(7, 7),
		dimnames=list (
			to=c('year', 'month', 'week', 'day',
			     'hour', 'minute', 'second'),
			from=c('year', 'month', 'week', 'day',
			       'hour', 'minute', 'second')))

	conv <- conversion[cbind (as.character(value),
				  as.character(unit(object)))]
	
	if (any (is.na(conv)))
		warning('some POSIXctp can not be converted due to incompatible units (year to hour for instance)')

	POSIXctp(duration(object) * ifelse(is.na(conv), 1, conv),
		 ifelse(is.na(conv), as.character (unit(object)),
				     as.character(value)))
} )

setMethod ('duration', 'POSIXctp', function(x, ...) x@duration)

format.POSIXctp <- function (x, ...) {
	ifelse( is.na(x@duration),
		'NA', 
		sprintf('%i %s%s',
			x@duration,
			as.character( x@unit ),
			ifelse(abs(x@duration)>1, 's', '')))
}

print.POSIXctp <- function(x, ...) print (format (x) )

setMethod ('show', 'POSIXctp', function(object) show (format (object) ))

tail.POSIXctp <- function (x, ...) tail(format(x, ...))

head.POSIXctp <- function (x, ...) head(format(x, ...))

summary.POSIXctp <- function (object, ...)
	summary(format(object, ...))

setMethod ('length', 'POSIXctp', function(x)length (x@duration))

'[<-.POSIXctp' <- function (x, i, value) {
	d <- duration (x)
	d[i] <- duration (value)
	u <- unit (x)
	u[i] <- unit (value)
	new ('POSIXctp', duration=d, unit=u)
}
'[.POSIXctp' <- function (x, i, ...) new ('POSIXctp', duration=duration (x)[i], unit=unit (x)[i])

as.POSIXctp.logical <- function (from, ...)
	if (is.na(from))
		POSIXctp(from) else
		stop ('Cannot coerce a logical to a POSIXctp.')

c.POSIXctp <- function(...){
	pers <- list(...)
	if (!all (sapply (pers, inherits, 'POSIXctp') ) )
		NextMethod('c')
	else
		new('POSIXctp',
		    duration=unlist(lapply(pers, duration)),
		    unit=POSIXt.units(unlist(lapply(pers, unit))))
}

Ops.POSIXctp <- function (e1, e2) {
	if (!inherits (e2, 'POSIXctp') ) return (NextMethod (.Generic) )
	if (!.Generic %in% c('==', '!=', '<=', '<', '>', '>='))
		return (NextMethod(.Generic) )

	l1 <- length(e1)
	l2 <- length(e2)
	if(l1 %% l2 != 0 & l2 %% l1 != 0) warning('longer object length is not',
						  ' a multipe of shorter',
						  ' object length')
	if (l1 > l2) e2 <- rep(e2, length.out=l1) else
	if (l1 < l2) e1 <- rep(e1, length.out=l2)
	
	# first we try to convert the POSIXctp with the bigger unit to
	# the unit of the other (if needed !)

	if(any( unit(e1) > unit(e2) )) suppressWarnings(
		unit(e1[unit(e1) > unit(e2)]) <- unit(e2)[unit(e1) > unit(e2)])

	if(any( unit(e2) > unit(e1) )) suppressWarnings(
		unit(e2[unit(e2) > unit(e1)]) <- unit(e1)[unit(e2) > unit(e1)])

	# then if two elements haven't the same unit, it means one is not
	# comparable to the other thus the result is NA. Otherwise, durations
	# are compared.

	# BUT if they are not comparable, they are not equal (or are different !!)

	if( .Generic == '==' ) {
		return( unit(e1)==unit(e2) & duration(e1)==duration(e2) )
	} else if( .Generic == '!=' ) {
		return( unit(e1)!=unit(e2) | duration(e1)!=duration(e2) )
	} else {
		return( ifelse(unit(e1) != unit(e2), NA, 
		       callGeneric(duration(e1), duration(e2))))
	}
}

setMethod('as.numeric', 'POSIXctp', function(x, ...) return( x@duration ))

setMethod('match', signature('POSIXctp', 'POSIXctp'),
	function(x, table, nomatch = NA_integer_, incomparables=NULL) {
		d <- match(x@duration, table@duration, NA, incomparables)
    		u <- match(unit(x), unit(table), NA, incomparables)

		ifelse(is.na(d) | is.na(u), nomatch,
		ifelse(d == u, d, nomatch))

	} )
setMethod('match', signature('POSIXctp', 'ANY'),
	function(x, table, nomatch = NA_integer_, incomparables=NULL)
	{
		if( length(unique(unit(x))) != 1 )
			stop("unit of table can't be guessed from x")
		table <- POSIXctp(table, unique(unit(x)))
		match(x, table, nomatch, incomparables)
	} )
setMethod('%in%', signature('POSIXctp', 'ANY'),
	function(x, table) match(x, table, 0L) > 0L)

unique.POSIXctp <- function(x, incomparables=FALSE, ...)
{
	u <- x[1]
	if(length(x) > 1)
	for( i in 2:length(x) ) {
		xi <- x[i]
		if(!any(duration(xi) == duration(u) &
			as.character(unit(xi)) == as.character(unit(u))))
			u <- c(u, xi)
	}
	return( u )
}

rep.POSIXctp <- function(x, ...) POSIXctp(rep(duration(x), ...), rep(unit(x), ...))
