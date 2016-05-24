# definition de la classe
#------------------------
setClass ('POSIXcti',
	representation (start='POSIXct', duration='integer'),
	validity=function(object) {
		if (	length (object@start) != length (object@duration) ) {
			return ("Slots 'start' and 'duration' must have same length.")
		} else return (TRUE)
	})

# constructeur
#-------------
POSIXcti <- function (start, end, timezone='UTC', ...)
{
	if (!inherits(start, 'POSIXct'))
		start <- as.POSIXct(start, tz=timezone)
	if (!inherits(end, 'POSIXct'))
		end <- as.POSIXct(end, tz=timezone)
	duration <- as.numeric(difftime(end, start, units='secs'))
	if (inherits(duration, 'numeric') & !inherits(duration, 'integer')) {
	#	warning('duration is not an integer. It will be coerced to')
		duration <- as.integer (duration)
	}
	return (new('POSIXcti', start=start, duration=duration))
}

# acces aux proprietes
#---------------------
setMethod ('duration', 'POSIXcti', function(x, ...) x@duration)

start.POSIXcti <- function (x, ...) {return (x@start) }
end.POSIXcti <- function (x, ...) start (x) + duration(x)

format.POSIXcti <- function (x, format='%Y-%m-%d %H:%M:%S', ...) {
	if (length (x) == 0) return ('POSIXcti()')
	tz <- attributes(x)$tzone
	s  <- start(x)
	e  <- end(x)
	tp <- paste(format(s, format=format),
		    format(e, format=format),
		    sep=' <-> ')
	return (tp)
}

print.POSIXcti <- function(x, ...) print (format (x) )
setMethod ('show', 'POSIXcti', function(object) show (format (object) ))

tail.POSIXcti <- function (x, ...) tail(format(x, ...))

head.POSIXcti <- function (x, ...) head(format(x, ...))

summary.POSIXcti <- function (object, ...)
	summary(format(object, ...))

setMethod ('length', 'POSIXcti', function(x)length (x@start))
'[<-.POSIXcti' <- function (x, i, value) {
	s <- start (x)
	s[i] <- start (value)
	d <- duration (x)
	d[i] <- duration (value)
	new ('POSIXcti', start=s, duration=d)
}
'[.POSIXcti' <- function (x, i, ...) new ('POSIXcti', start=start(x)[i], duration=duration (x)[i])

as.POSIXcti.logical <- function (from, ...)
	if (is.na(from))
		new ('POSIXcti', start=as.POSIXct(NA), duration=as.integer (NA)) else
		stop ('Cannot coerce a logical to a POSIXcti.')

'%intersect%.POSIXcti' <- function(i1, i2) {
	s1 <- start (i1)#attr (i1, 'start')
	s2 <- start (i2)#attr (i2, 'start')
	e1 <- end(i1)#i1
	e2 <- end(i2)#i2
	
	s <- mapply (max, s1, s2)
	e <- mapply (min, e1, e2)

	if (any (s > e) )
		warning ("Some intersections are empty. NA's introduced.")

	interval <- POSIXcti(as.POSIXct(s, origin=origin), as.POSIXct(e, origin=origin) )
	interval[s > e] <- as.POSIXcti(NA)
	return (interval)
}

'%included%.POSIXcti' <- function(i1, i2) {
	mapply ('&', SIMPLIFY=TRUE,
		mapply ('>=', SIMPLIFY=FALSE, start(i1), start(i2)),
		mapply ('<=', SIMPLIFY=FALSE, end(i1), end(i2))
		)
}


# '%union%' <- function(i1, i2) {}

c.POSIXcti <- function(...){
	pers <- list(...)
	if (!all (sapply (pers, inherits, 'POSIXcti') ) ) {
		NextMethod('c')
	} else {
		if( length(unique(sapply(pers, function(x) attributes(start(x))$tz[1]))) > 1 )
			warning( "timezone are differents. The first is used." )
		s <- do.call(c, lapply(pers, start))
		new('POSIXcti', start=s, duration=unlist(lapply(pers, duration)))
	}
}

Ops.POSIXcti <- function (e1, e2) {
	if (!inherits (e2, 'POSIXcti') ) return (NextMethod (.Generic) )
	if (!.Generic %in% c('==', '!=', '<=', '<', '>', '>='))# stop (sprintf ("%s not implemented for 'period' objects"), .Generic)
		#                 stop(gettextf("'%s' not defined for \"period\" objects", 
		#                                           .Generic), domain = NA)
		return (NextMethod(.Generic) )
	if (.Generic == '==')
		return (start(e1) == start(e2) & duration(e1) == duration(e2))
	if (.Generic == '!=')
		return (!e1 == e2)
	if (.Generic == '<') {
		return (end(e1) <= start(e2))
	}
	if (.Generic == '<=') {
		return (start(e1) <= start(e2) & end(e1) <= end(e2))
	}
	if (.Generic == '>')
		return (e2 < e1)
	if (.Generic == '>=')
		return (e2 <= e1)
}

setMethod('match', signature('POSIXcti', 'POSIXcti'),
	function(x, table, nomatch = NA_integer_, incomparables=NULL)
	{
		m <- match(x@start, table@start, nomatch, incomparables)
		m[table@duration[m]!=x@duration] <- nomatch
		m
	} )
setMethod('%in%', signature('POSIXcti', 'POSIXcti'),
	function(x, table) match(x, table, 0L) > 0L)

unique.POSIXcti <- function(x, incomparables=FALSE, ...)
{
	u <- x[1]
	if(length(x) > 1)
	for( i in 2:length(x) ) {
		xi <- x[i]
		if(!any(duration(xi) == duration(u) &
			start(xi) == start(u)))
			suppressWarnings( u <- c(u, xi) )
	}
	return( u )
}

rep.POSIXcti <- function(x, ...)
	new('POSIXcti', start=rep(start(x), ...), duration=rep(duration(x), ...))
