# definition de la classe
#------------------------
setClass (Class = 'POSIXst', 
	  representation = representation (
		subtime='integer', unit='factor', of='factor',
	       	timezone='character'),
	  prototype = prototype (
		subtime=integer(), unit=POSIXt.units('second'),
		of=POSIXt.units('hour'), timezone='UTC'),
	  validity=function(object) {
		if( length(object@unit) != 1 )
			stop("length of 'unit' should be equal to 1")
		if( length(object@of) != 1 )
			stop("length of 'of' should be equal to 1")
		if( !object@unit %in% POSIXt.units() )
			  stop(sprintf("'unit' should one of %s",
			  paste(as.character(POSIXt.units()), collapse=", ")))
		if( !object@of %in% POSIXt.units() )
			  stop(sprintf("'of' should one of %s",
			  paste(as.character(POSIXt.units()), collapse=", ")))
		if( object@unit >= object@of )
			stop("'unit' should be a subdivision time of 'of'")

		uc <- as.character( object@unit )
		oc <- as.character( object@of )

		if( uc == 'week' & oc == 'month' )
			stop( "week of month not implemented")
		if( uc == 'week' & oc == 'year' )
			stop( "week of year not implemented")

		if( uc=='year' ) {
			val.min <- -Inf
			val.max <- +Inf
		} else if( uc=='month' ) {
			val.min <- 0
			val.max <- 11
		} else if( uc=='day' & oc=='month' ) {
			val.min <- 1
			val.max <- 31
		} else {
			val.min <- 0

			bases <- switch(oc, year=c(60, 60, 24, 366),
					    month=c(60, 60, 24, 31),
					    week=c(60, 60, 24, 7),
					    day=c(60, 60, 24),
					    hour=c(60, 60),
					    minute=c(60))
			sb <- switch(uc, second=1, minute=2, hour=3, 4)
			val.max <- prod(bases[sb:length(bases)])-1

			if( uc=='second' ) val.max <- val.max+2
		}
		if( is.null(oc) ) oc <- ''
		if( !all(val.min <= object@subtime & object@subtime <= val.max ) )
			stop(sprintf('For %s%s, object@subtime should be between %i and %i',
				uc, ifelse(oc=='', '', sprintf(' of %s', oc)), val.min, val.max))

		  return (TRUE)
	  } )

#--------------
# constructors
# see the POSIXst file in which all constructors are
# defined

#---------------------
# properties accessors

unit.POSIXst <- function(x, ...) x@unit

of.POSIXst <- function(x, ...) x@of

timezone.POSIXst <- function(object) object@timezone

#-----------------
# printing methods

format.POSIXst <- function (x, format=NULL, ...)
{
	if( is.null(format) )
	{
		if( x@unit == POSIXt.units('month') )
			format <- '%B' else
		if( x@unit == POSIXt.units('day') &
		    x@of == POSIXt.units('week') )
			format <- '%A' else
		if( x@unit == POSIXt.units('year') )
			format <- '%s %v' else
			format <- '%v%p %s of %m'
	}

	result <- rep(format, length(x))

	# case of days of week
	if( any( sapply(c('%a', '%A'), grepl, format) ) ) 
		if( x@unit != POSIXt.units('day') &&
		    x@of != POSIXt.units('week') ) {
			stop('%a, %A format can only be used with month.')
		} else {
			tmp <- as.POSIXct(sprintf('2011-12-%02i', 25:31))[
								x@subtime+1]
			result <- mapply(
				function(result, tmp)
				{
					gsub('%a', format(tmp, '%a'),
					gsub('%A', format(tmp, '%A'),
					result))
				},
				result, tmp, USE.NAMES=FALSE)
		}

	# case of month
	if( any( sapply(c('%b', '%B'), grepl, format) ) )
		if( x@unit != POSIXt.units('month') ) {
			stop( '%b, %B format can only be used with month.' )
		} else {
			tmp <- as.POSIXct(sprintf('1970-%02i-01', x@subtime+1))
			result <- mapply(function(result, tmp) {
				 gsub('%b', format(tmp, '%b'),
				 gsub('%B', format(tmp, '%B'),
				      result))
				},
				result, tmp, USE.NAMES=FALSE)
		}

	result <- mapply(function(result, p, v, s, m, t)
	 		 gsub('%p', p,
	 		 gsub('%v', v,
	     		 gsub('%s', s,
			 gsub('%m', m,
			 gsub('%t', t, result))))),
			 result,
			sapply(sprintf('a%i', x@subtime%%10), switch,
			       a1='st', a2='nd', a3='rd', 'th'),# p
			 x@subtime, 				# v
			 as.character(x@unit),			# s
			 as.character(x@of),			# m
			 x@timezone[1],				# t
			 USE.NAMES=FALSE)
	result
}

print.POSIXst <- function (x, ...)
	print(format(x, ...))

setMethod ('show', 'POSIXst', function (object) print (object))

tail.POSIXst <- function (x, ...) tail(format(x, ...))

head.POSIXst <- function (x, ...) head(format(x, ...))

summary.POSIXst <- function (object, ...) summary(format(object, ...))

#-----------------
# common accessors
'[.POSIXst' <- function(x, i)
{
	if(missing(i)) i <- seq_len(length(x))
	x@subtime <- x@subtime[i]
	validObject(x)
	return(x)
}

'[<-.POSIXst' <- function(x, i, value)
{
	if( !inherits(value, 'POSIXst') )
		stop("'value' is not a POSIXst object")
	if( unit(x) != unit(value) )
		stop("'x' and 'value' must have same unit")
	if( of(x) != of(value) )
		stop("'x' and 'value' must have same of")
	if( timezone(x) != timezone(value) )
		stop("'x' and 'value' must have same timezone")
	if(missing(i)) i <- seq_len(length(x))
	if( length(i) != length( value ) )
		stop("'i' and 'value' must have same length.")
	x@subtime[i] <- value
	validObject(x)
	return(x)
}
setMethod (f='length', signature='POSIXst',
	   definition=function(x) length (x@subtime))

# Math

# manipulation
#-------------
c.POSIXst <- function(...)
{
	pers <- list(...)
	if (!all (sapply (pers, inherits, 'POSIXst') ) )
		NextMethod('c') else {
		us <- unique(lapply( pers, unit ))
		if( length(us) != 1)
			stop("all POSIXst must have same 'unit' to be used with the 'c' function")

		os <- unique(lapply( pers, of ))
		if( length(os) != 1)
			stop("all POSIXst must have same 'of' to be used with the 'c' function")

		tzs <- unique(lapply( pers, function(x) timezone(x)[1]))
		if( length(tzs) != 1)
			stop("all POSIXst must have same 'timezone' to be used with the 'c' function")

		new('POSIXst', subtime=unlist(lapply(pers, slot, 'subtime')),
	   unit=us[[1]], of=os[[1]], timezone=tzs[[1]])
	}
}
Ops.POSIXst <- function (e1, e2) {
	if( !inherits (e2, 'POSIXst') ) return (NextMethod (.Generic) )
	if( !.Generic %in% c('==', '!=', '<=', '<', '>', '>=') )
		return (NextMethod(.Generic) )
	if( .Generic == '==' )
	{
		return( unit(e1) == unit(e2) & of(e1) == of(e2) & timezone(e1) == timezone(e2) &
			e1@subtime == e2@subtime )
	}
	if( .Generic == '!=' )
		return (!e1 == e2)

	if( unit(e1) != unit(e2) )
		stop( "e1 and e2 must have the same 'unit' to be compared" )
	if( of(e1) != of(e2) )
		stop( "e1 and e2 must have the same 'of' to be compared" )
	if( timezone(e1) != timezone(e2) )
		stop( "e1 and e2 must have the same 'timezone' to be compared" )

	do.call( .Generic, list(e1=e1@subtime, e2=e2@subtime) )
}

setMethod('as.numeric', 'POSIXst', function(x, ...) return( x@subtime ))

setMethod('match', signature('POSIXst', 'POSIXst'),
	function(x, table, nomatch = NA_integer_, incomparables=NULL)
	{
		if( unit(x) != unit(table) ) return(rep(nomatch, length(x)))
		if( of(x) != of(table) ) return(rep(nomatch, length(x)))
		match(x@subtime, table@subtime, nomatch, incomparables)
	} )
setMethod('match', signature('POSIXst', 'ANY'),
	function(x, table, nomatch = NA_integer_, incomparables=NULL)
	{
		table <- POSIXst(table, unit(x), of(x))
		match(x, table, nomatch, incomparables)
	} )
setMethod('%in%', signature('POSIXst', 'ANY'),
	function(x, table) match(x, table, 0L) > 0L)

unique.POSIXst <- function(x, incomparables=FALSE, ...)
	x[!duplicated(x@subtime)]

duplicated.POSIXst <- function(x, incomparables=FALSE, ...)
	duplicated(x@subtime)

rep.POSIXst <- function(x, ...)
	POSIXst(rep(x@subtime, ...), unit(x), of(x), timezone(x)[1])

seq.POSIXst <- function(from, to, ...) {
	if( unit(from) != unit(to) | of(from) != of(to) )
		stop("'from' and 'to' must have same units")
	if( timezone(from) != timezone(to) )
		warning("'timezone' of 'from' and 'to' are different, the timezone of 'from' is used")

	return(POSIXst(seq(as.numeric(from), as.numeric(to), ...),
		       unit(from), of(from), timezone(from)))
}
