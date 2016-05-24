# definition de la classe
#------------------------
setClass (Class = 'TimeInstantDataFrame', 
	  representation = representation (instant='POSIXct', timezone='character', data='data.frame'),
	  prototype = prototype (instant=as.POSIXct(character(), timezone='UTC'),
				 data=data.frame()),
	  validity=function(object) {
		  if (length (when (object)) != nrow (object))
			  stop ("In a 'TimeInstantDataFrame, 'data' must have a number of rows as long as 'instant'.")
		  return (TRUE)
	  })

# constructeurs
#--------------
TimeInstantDataFrame <- function (when, timezone='UTC', data=NULL, ...) {
	if (is.character (when) ) when <- as.POSIXct (when, timezone)
	if (is.null (data)) data <- data.frame (matrix (NA, ncol=0, nrow=length(when) ) )
	new ('TimeInstantDataFrame', instant=when, timezone=timezone, data=data)
}

# Create a regular TimeInstantDataFrame from scratch
RegularTimeInstantDataFrame <- function (from, to, by, timezone='UTC', data=NULL) {
	if (is.character (from) ) from <- as.POSIXct (from, timezone)
	if (is.character (by) ) by <- POSIXctp(unit=by)
	if (missing (to))
		to <- from + (nrow(data) - 1) * by
	if (is.character (to) ) to <- as.POSIXct (to, timezone)
	if (!inherits (by, 'POSIXctp') )
		stop ("'by' should be coercible to a 'POSIXctp'.")

	if (as.character(unit(by)) == 'year') {
		nb <- year(to) - year(from) + 
			ifelse(second(to, of='year') == 0, 0, 1)
	} else if (as.character(unit(by)) == 'month') {
		nb <- (year(to) - year(from))*12 + as.numeric(month(to)) - as.numeric(month(from)) + 
			ifelse(second(to, of='month') == 0, 0, 1)
	} else {
		u <- switch (as.character(unit(by)), second='secs', minute='mins',
						     hour='hours', day='days')
		nb <- as.numeric (difftime(to, from, units=u))
		nb <- ceiling (nb/duration(by))
	}
	when <- from + 0:nb * by
	tk <- !is.na(when) & (when >= from & when <= to) 
	when <- when[tk]

	if (is.null (data)) data <- data.frame (matrix (NA, ncol=0, nrow=length(when) ) )
	new ('TimeInstantDataFrame', instant=when, timezone=timezone, data=data)
}

# definition des accesseurs de l'objet
#-------------------------------------

setMethod (f='when', signature='TimeInstantDataFrame',
definition=function(x, ...)
	return(as.POSIXct(as.POSIXlt(x@instant, timezone(x)))) )

setMethod (f='timezone', signature='TimeInstantDataFrame',
	   definition=function(object) return(object@timezone[1]) )

setMethod (f='timezone<-', signature='TimeInstantDataFrame',
definition=function(object, value) {
	object@timezone <- value
	object@instant  <- as.POSIXct(as.POSIXlt( object@instant, value ))
	return(object) } )

# mise en forme pour / et affichage
#----------------------------------
print.TimeInstantDataFrame <- function (x, tz=NULL, ...) {
	if (is.null (tz) ) tz <- timezone(x)
	print(data.frame (when=format (when(x), tz=tz, usetz=TRUE), x@data) )
}
setMethod ('show', 'TimeInstantDataFrame',
	   function (object) print (object, timezone(object))
)
		   #                    print(data.frame (when=when(object), object@data) ), tz=timezone(object))
tail.TimeInstantDataFrame <- function (x, tz, ...) {
	if (missing (tz) ) tz <- x@timezone
	tail(data.frame (when=format (when(x), tz=tz, usetz=TRUE), x@data), ...)
}
head.TimeInstantDataFrame <- function (x, tz, ...) {
	if (missing (tz) ) tz <- x@timezone
	head(data.frame (when=format (when(x), tz=tz, usetz=TRUE), x@data), ...)
}
summary.TimeInstantDataFrame <- function (object, ...)
		summary(data.frame (when=when(object), object@data), ...)
# format

# defintion des accesseurs aux donnees
#-------------------------------------
'[.TimeInstantDataFrame' <- function(x, i, j, drop=FALSE) {
	n.args <- nargs() - hasArg(drop)
	if (missing (j) & n.args==2) {
		j <- i
		i <- rep(TRUE, nrow(x))
	}

	if(missing(i)) i <- rep(TRUE, nrow(x))

	# for i = 'YYYY-MM-DD HH:MM:SS tz'
	if (!missing(i) && length(i) == 1 && (is.character(i) || inherits(i, 'POSIXt'))) {
		if (is.character(i)) {
			di <- strsplit(i, ' ')[[1]]
			if (length(di) == 2 && !grepl('..:..:..', di[2])) {
				di[3] <- di[2]
				di[2] <- ''
			}
			if(is.na(di[3])) di[3] <- timezone(x)

			di <- try(as.POSIXct(paste(di[1], di[2]), di[3]), TRUE)
		} else di <- i

		if ( !inherits(di, 'try-error') ) i <- when(x) >= di
	}
	# for j = 'YYYY-MM-DD HH:MM:SS tz'
	if (!missing(j) && length(j) == 1 && (is.character(j) || inherits(j, 'POSIXt'))) {
		if (is.character(j)) {
			di <- strsplit(j, ' ')[[1]]
			if (length(di) == 2 && !grepl('..:..:..', di[2])) {
				di[3] <- di[2]
				di[2] <- ''
			}
			if(is.na(di[3])) di[3] <- timezone(x)

			di <- try(as.POSIXct(paste(di[1], di[2]), di[3]), TRUE)
		} else di <- j

		if ( !inherits(di, 'try-error') ){
			i <- i & when(x) <= di
			j <- names(x)
		}
	}

	y <- new ('TimeInstantDataFrame', 
	     instant = when(x)[i, drop=drop],
	     data = x@data[i, j, drop=drop],
	     timezone=timezone(x)[1])
	validObject(y)
	return(y)
}
setMethod (f='[[', signature='TimeInstantDataFrame',
	   definition=function(x, i, ...) {
		   '[[.data.frame'(x@data, i, ...)
	   })
setMethod (f='$', signature='TimeInstantDataFrame',
	   definition=function(x, name) {
		   do.call ('$', list(x=x@data, name=name))
	   })

'[<-.TimeInstantDataFrame' <- function(x, i, j, value) {
	n.args <- nargs()
	if (missing (j) & n.args==3) {
		j <- i
		i <- seq_len(nrow(x))
	}
	if(missing(i)) i <- seq_len(nrow(x))
	x@data[i,j] <- value
	validObject(x)
	return(x)
}
'[[<-.TimeInstantDataFrame' <- function(x, i, j, value) {
   if (missing (j) )
	   x@data[[i]] <- value else
	   x@data[[i,j]] <- value
   validObject(x)
   return(x)
}
setMethod (f='$<-', signature='TimeInstantDataFrame',
	   definition=function(x, name, value) {
		   x@data <- "$<-.data.frame"(x@data, name, value)
		   validObject(x)
		   return(x)
	   })

setMethod (f='dim', signature='TimeInstantDataFrame',
	   definition=function(x) dim (x@data))
setMethod (f='length', signature='TimeInstantDataFrame',
	   definition=function(x) length (x@data))
setMethod (f='nrow', signature='TimeInstantDataFrame',
	   definition=function(x) nrow (x@data))
setMethod (f='ncol', signature='TimeInstantDataFrame',
	   definition=function(x) ncol (x@data))
row.names.TimeInstantDataFrame <- function(x) row.names (x@data)
'row.names<-.TimeInstantDataFrame' <- function(x, value) {
		   row.names (x@data) <- value
		   x
	   }
setMethod (f='names', signature='TimeInstantDataFrame',
	   definition=function(x) names (x@data))
setMethod (f='names<-', signature='TimeInstantDataFrame',
	   definition=function(x, value) {
		   names (x@data) <- value
		   x
	   } )

# Math

# manipulation
#-------------
# fonction réalisée en S3 pour ne pas imposer de 'signature'
rbind.TimeInstantDataFrame <- function (...)
{
	dots <- list (...)
	names(dots) <- NULL
	if (!all (sapply (dots, inherits, 'TimeInstantDataFrame')))
		stop ("all arguments must be 'TimeInstantDataFrame'")
	instant <- as.POSIXct (unlist (lapply (dots, when) ),
			       origin=timetools::origin)
	df <- do.call("rbind", lapply(dots, function(x) x@data) )
	tz <- timezone (dots[[1]])
	if (!all (tz == sapply (dots, timezone)))
		warning ("Not all timezone are identical. Timezone of the first object is used.")
	new('TimeInstantDataFrame', instant=instant,
	    timezone=timezone (dots[[1]]), data=df)
}
# cbind # a faire eventuellement entre un Time*DataFrame et une data.frame
merge.TimeInstantDataFrame <- function(x, y, by, all=TRUE, tz='UTC', ...)
{
	instant.vec <- list (when(x), when(y))
	x.data <- data.frame (instant=format (when(x),
					      format='%Y-%m-%d %H:%M:%S',
					      tz='UTC'),
			      x@data)
	y.data <- data.frame (instant=format (when(y),
					      format='%Y-%m-%d %H:%M:%S',
					      tz='UTC'),
			      y@data)
	z <- merge (x.data, y.data, by=unique (c('instant', by) ), all=all, ...)
	z <- new ('TimeInstantDataFrame',
     		  instant=as.POSIXct(z$instant, tz='UTC'),
     		  data=z[setdiff(names(z), c('instant'))],
     		  timezone=tz)
	timezone(z) <- tz
	return (z)
}

setMethod ('lapply', signature('TimeInstantDataFrame', 'ANY'),
	   function (X, FUN, ...)
	   {
		   res <- lapply (data.frame(X), FUN, ...)
		   if (all (sapply (res, length) == nrow(X))) {
			   X@data <- data.frame (res[names(X)])
		   } else if (all (sapply (res, length) == 1)) {
			   X <- new ('TimeIntervalDataFrame',
				     start=min(when(X)), end=max(when(X)),
				     timezone=timezone(X),
				     data=data.frame (res))
		   } else {
			   stop ("try to apply inadequate function over TimeInstantDataFrame.")
		   }
		   return (X)
	   } )

# acces/modification de certaines propriétés
#-------------------------------------------
setMethod (f='regular', signature='TimeInstantDataFrame',
definition=function(x, ...) {
	len <- length(unique(difftime(when(x)[-1], when(x)[-nrow(x)])))
	return(length(len) == 1)
} )

# transformateur de classe
#-------------------------
setAs ('TimeInstantDataFrame', 'data.frame',
       function(from) data.frame (instant=when(from), from@data) )

as.data.frame.TimeInstantDataFrame <- function (x, row.names=NULL, optional=FALSE, include.dates=FALSE, ...) {
	if (include.dates)
		return (data.frame (date=when (x), x@data) ) else
		return (x@data)
}

as.TimeIntervalDataFrame.TimeInstantDataFrame <- function(from, period, ...) {
	if (nrow(from) == 0) {
		to <- TimeIntervalDataFrame (
			   start=character(0), end=character(0), 
			   timezone=timezone(from), data=from@data)

	} else {
		if (missing(period)) {
			if (regular(from)) {
				period <- as.numeric(difftime (when(from)[2], when(from)[1],
											   units='secs'))
			} else {
				stop ("'period' must be of class 'period' or 'from' should be ",
					  "at least 'regular'.")
			}
			period <- POSIXctp (period, 'second')
		}

		to <- new ('TimeIntervalDataFrame',
				   start=when(from), end=when(from)+period, 
				   timezone=timezone(from), data=from@data)
	} # fin du if sur nrow(from) == 0

	validObject(to)
	return (to)
}

as.SubtimeDataFrame.TimeInstantDataFrame <- function(x, unit, of, FUN=NULL, ...)
{
	st <- POSIXst(x, unit, of)
	to <- data.frame( x )
	if( !is.null(FUN) )
	{
		u <- unit(st)
		o <- of(st)
		tz <- timezone(st)
		st <- as.numeric(format( st, "%v" ))
		to <- split (to, st)
		to <- lapply (to, sapply, FUN, ...)
		st <- POSIXst(as.numeric(names(to)), u, o, tz)
		to <- t(data.frame (to))
		rownames (to) <- NULL
	}

	to <- new ('SubtimeDataFrame', when=st, data=data.frame (to))
	validObject(to)
	return (to)
}
