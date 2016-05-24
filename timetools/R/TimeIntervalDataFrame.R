# definition de la classe
#------------------------
setClass (Class = 'TimeIntervalDataFrame', 
	  representation = representation (start='POSIXct', end='POSIXct',
					   timezone='character', data='data.frame'),
	  prototype = prototype (start=as.POSIXct(character(), tz='UTC'),
				 end=as.POSIXct(character(), tz='UTC'),
				 timezone='UTC',
				 data=data.frame()),
	  validity=function(object) {
		  if (length (object@start) != length (object@end))
			  stop ("'start' and 'end' must have the same length in a 'TimeIntervalDataFrame'.")
		  if (length (object@end) != nrow (object@data))
			  stop ("In a 'TimeIntervalDataFrame, 'data' must have a number of rows as long as 'start' and 'end'.")
		  return (TRUE)
	  } )

# constructeurs
#--------------
TimeIntervalDataFrame <- function (start, end=NULL, timezone='UTC', data=NULL, period=NULL, ...) {
	# cas avec period

	if(!is.null(period)) {
		if(length(start) != 1 & length(end) != 1)
			stop("both 'start' and 'end' arguments must have a length of 1.")
			
		if (length (period) > 1) {
			warning ('Only the first given period is used as \'to\'.')
			period <- period[1]
		}

		if( is.character(start) ) start <- as.POSIXct (start, timezone)
		if( is.character(end) ) end <- as.POSIXct (end, timezone)
		if( is.character(period) ) period <- POSIXctp( period )
		tz <- if(is.null(timezone)) attributes(start)$tzone[1] else timezone

		# construction de ls structure qui va servir 

		u <- as.character (unit(period))
		if (u == 'second') {
			s <- trunc (start, 'secs')
		} else if (u == 'minute') {
			s <- trunc (start, 'mins')
		} else if (u == 'hour') {
			s <- trunc (start, 'hours')
		} else if (u == 'day') {
			s <- trunc (start, 'days')
		} else if (u == 'month') {
			s <- as.POSIXct(
				sprintf('%s-01', format(start, '%Y-%m')),
				tz)
		} else if (u == 'year') {
			s <- as.POSIXct(
				sprintf('%s-01-01', format(start, '%Y')),
				tz)
		}
		
		s <- as.POSIXct(s)
		
		if (u == 'year') {
			e <- end
			nb <- year(e) - year(s) +
			ifelse(second(e, of='year') == 0, 0, 1)
		} else if (u == 'month') {
			e <- end
			nb <- (year(e) - year(s))*12 +
				as.numeric(month(e)) -
				as.numeric(month(s)) +
			ifelse(second(e, of='month') == 0, 0, 1)
		} else {
			u <- switch (u, second='secs', minute='mins',
				     hour='hours', day='days')
			nb <- as.numeric (difftime(end, s, units=u))
			nb <- ceiling (nb/duration(period))
		}
		
		e <- s+as.numeric(nb) * period
		
		result <- RegularTimeIntervalDataFrame(
			s, e, by=period, timezone=tz)
		return( result )
	}

	# cas classique

	if (is.null (end) ) {
		end <- start[-1]
		start <- start[-length(start)]
	}
	if (is.character (start) ) start <- as.POSIXct (start, timezone)
	if (is.character (end) ) end <- as.POSIXct (end, timezone)
	if (is.null (data)) data <- data.frame (matrix (NA, ncol=0, nrow=length(start) ) )
	new ('TimeIntervalDataFrame', start=start, end=end,
	     timezone=timezone, data=data)
}

# Create a regular TimeIntervalDataFrame from scratch
RegularTimeIntervalDataFrame <- function (from, to, by, period, timezone='UTC', data=NULL) {
	if (is.character (from) ) from <- as.POSIXct (from, timezone)
	if (is.character (by) ) by <- POSIXctp(unit=by)
	if (missing (to))
		to <- from + nrow(data) * by
	if (is.character (to) ) to <- as.POSIXct (to, timezone)
	if (!inherits (by, 'POSIXctp') )
		stop ("'by' should be coercible to a 'POSIXctp'.")
	if (missing (period) ) {
		period <- by
	} else if (is.character (period) ) {
		period <- POSIXctp(unit=period)
	}
	if (!inherits (period, 'POSIXctp') )
		stop ("'period' should be coercible to a 'POSIXctp'.")

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
	}
	nb <- ceiling (nb/duration(by))
	start <- from + 0:(nb-1) * by
	end <- start + period
	tk <- !is.na(start) & !is.na(end) &
		((start >= from & start <= to) | (start >= from & end <= to))
	start <- start[tk]
	end <- end[tk]

	if (is.null (data)) data <- data.frame (matrix (NA, ncol=0, nrow=length(start) ) )
	new ('TimeIntervalDataFrame', start=start, end=end,
	     timezone=timezone, data=data)
}

# definition des accesseurs de l'objet
#-------------------------------------

start.TimeIntervalDataFrame <- function(x, ...) return(as.POSIXct(as.POSIXlt(x@start, timezone(x))))
end.TimeIntervalDataFrame <- function(x, ...) return(as.POSIXct(as.POSIXlt(x@end, timezone(x))))

setMethod (f='timezone', signature='TimeIntervalDataFrame',
	   definition=function(object) return(object@timezone[1]) )

setMethod (f='timezone<-', signature='TimeIntervalDataFrame',
definition=function(object, value) {
	start <- as.POSIXct(as.POSIXlt( object@start, value ))
	end   <- as.POSIXct(as.POSIXlt( object@end, value ))

	new ('TimeIntervalDataFrame',
	     start=start, end=end,
	     timezone=value, data=data.frame (object)
	     ) } )

setMethod (f='interval', signature='TimeIntervalDataFrame',
definition=function(x, ...) return(POSIXcti (start(x), end(x)) ) )

setMethod (f='when', signature='TimeIntervalDataFrame',
definition=function(x, ...) return(POSIXcti (start(x), end(x)) ) )

setMethod (f='period', signature='TimeIntervalDataFrame',
definition=function(x, ...) {
	if (!homogeneous(x))
		stop ("x should be homogeneous to have a 'period'")

	if (!continuous(x))
		stop ("x should be continuous to have a 'period'")
	
	res <- difftime(end(x), start(x), units='secs')
	res <- POSIXctp (unique(res), unit=rep ('second', length(unique(res))))
	return (res)
} )

# mise en forme pour / et affichage
#----------------------------------
print.TimeIntervalDataFrame <- function (x, tz=NULL, ...) {
	if (is.null (tz) ) tz <- timezone(x)
	print(data.frame (start=format (start(x), tz=tz, usetz=TRUE),
			 end=format (end(x), tz=tz, usetz=TRUE),
			 x@data), ...)
}
setMethod ('show', 'TimeIntervalDataFrame',
	   function (object)
		   print (object, timezone(object)) )
tail.TimeIntervalDataFrame <- function (x, tz, ...) {
	if (missing (tz) ) tz <- x@timezone
	tail(data.frame (start=format (start(x), tz=tz, usetz=TRUE),
			 end=format (end(x), tz=tz, usetz=TRUE),
			 x@data), ...)
}
head.TimeIntervalDataFrame <- function (x, tz, ...) {
	if (missing (tz) ) tz <- x@timezone
	head(data.frame (start=format (start(x), tz=tz, usetz=TRUE),
			 end=format (end(x), tz=tz, usetz=TRUE),
			 x@data), ...)
}
summary.TimeIntervalDataFrame <- function (object, ...)
	summary(data.frame (start=start(object), end=end(object), object@data), ...)
# format

# defintion des accesseurs aux donnees
#-------------------------------------
'[.TimeIntervalDataFrame' <- function(x, i, j, drop=FALSE) {
	n.args <- nargs() - hasArg(drop)
	# call : x[i]
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

		if ( !inherits(di, 'try-error') ) i <- start(x) >= di
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
			i <- i & end(x) <= di
			j <- names(x)
		}
	}

	y <- new ('TimeIntervalDataFrame', 
	     start = start(x)[i, drop=drop],
	     end = end(x)[i, drop=drop],
	     timezone = x@timezone,
	     data = x@data[i, j, drop=drop])
	validObject(y)
	return(y)
}
setMethod (f='[[', signature='TimeIntervalDataFrame',
	   definition=function(x, i, ...) {
		   '[[.data.frame'(x@data, i, ...)
	   })
setMethod (f='$', signature='TimeIntervalDataFrame',
	   definition=function(x, name) {
		   do.call ('$', list(x=x@data, name=name))
	   })

'[<-.TimeIntervalDataFrame' <- function(x, i, j, value) {
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
'[[<-.TimeIntervalDataFrame' <- function(x, i, j, value) {
   if (missing (j) )
	   x@data[[i]] <- value else
	   x@data[[i,j]] <- value
   validObject(x)
   return(x)
}
setMethod (f='$<-', signature='TimeIntervalDataFrame',
	   definition=function(x, name, value) {
		   x@data <- "$<-.data.frame"(x@data, name, value)
		   validObject(x)
		   return(x)
	   })

setMethod (f='dim', signature='TimeIntervalDataFrame',
	   definition=function(x) dim (x@data))
setMethod (f='length', signature='TimeIntervalDataFrame',
	   definition=function(x) length (x@data))
setMethod (f='nrow', signature='TimeIntervalDataFrame',
	   definition=function(x) nrow (x@data))
setMethod (f='ncol', signature='TimeIntervalDataFrame',
	   definition=function(x) ncol (x@data))
row.names.TimeIntervalDataFrame <- function(x) row.names (x@data)
'row.names<-.TimeIntervalDataFrame' <- function(x, value) {
		   row.names (x@data) <- value
		   x
	   }
setMethod (f='names', signature='TimeIntervalDataFrame',
	   definition=function(x) names (x@data))
setMethod (f='names<-', signature='TimeIntervalDataFrame',
	   definition=function(x, value) {
		   names (x@data) <- value
		   x
	   } )

# Math

# manipulation
#-------------
# fonction réalisée en S3 pour ne pas imposer de 'signature'
rbind.TimeIntervalDataFrame <- function (...)
{
	dots <- list (...)
	names(dots) <- NULL
	if (!all (sapply (dots, inherits, 'TimeIntervalDataFrame')))
		stop ("all arguments must be 'TimeIntervalDataFrame'")
	start <- as.POSIXct (unlist (lapply (dots, start) ), origin=timetools::origin)
	end <- as.POSIXct (unlist (lapply (dots, end) ), origin=timetools::origin)
	df <- do.call("rbind", lapply(dots, function(x) x@data) )
	tz <- timezone (dots[[1]])
	if (!all (tz == sapply (dots, timezone)))
		warning ("Not all timezone are identical. Timezone of the first object is used.")
	new('TimeIntervalDataFrame', start=start, end=end,
	     timezone=tz, data=df)
}
# cbind # a faire eventuellement entre un Time*DataFrame et une data.frame
merge.TimeIntervalDataFrame <- function(x, y, by, all=TRUE, tz='UTC', ...) {
	if (!inherits(y, 'TimeIntervalDataFrame'))
		stop ("'y' must be a 'TimeIntervalDataFrame'.")
	if (missing (by) ) by <- intersect (names (x), names(y))
	start.vec <- list (start(x), start(y))
	end.vec <- list (end(x), end(y))
	x.data <- data.frame (start=format (start(x),
					    format='%Y-%m-%d %H:%M:%S',
					    tz='UTC'),
			      end=format (end(x),
					  format='%Y-%m-%d %H:%M:%S',
					  tz='UTC'),
			      x@data)
	y.data <- data.frame (start=format (start(y),
					    format='%Y-%m-%d %H:%M:%S',
					    tz='UTC'),
			      end=format (end(y),
					  format='%Y-%m-%d %H:%M:%S',
					  tz='UTC'),
			      y@data)
	z <- merge (x.data, y.data, by=unique (c('start', 'end', by) ), all=all, ...)
	z <- new ('TimeIntervalDataFrame',
	     start=as.POSIXct(z$start, tz='UTC'),
	     end=as.POSIXct(z$end, tz='UTC'),
	     timezone='UTC',
	     data=z[setdiff(names(z), c('start', 'end'))])
	timezone(z) <- tz
	return (z)
}

setMethod ('lapply', signature('TimeIntervalDataFrame', 'ANY'),
	   function (X, FUN, ...)
	   {
		   res <- lapply (data.frame(X), FUN, ...)
		   if (all (sapply (res, length) == nrow(X))) {
			   X@data <- data.frame (res[names(X)])
		   } else if (all (sapply (res, length) == 1)) {
			   X <- new ('TimeIntervalDataFrame',
				     start=min(start(X)), end=max(end(X)),
				     timezone=timezone(X),
				     data=data.frame (res))
		   } else {
			   stop ("try to apply inadequate function over TimeIntervalDataFrame.")
		   }
		   return (X)
	   } )
# acces/modification de certaines propriétés
#-------------------------------------------
# tous les intervalles sont de même durée
setMethod (f='homogeneous', signature='TimeIntervalDataFrame',
definition=function(x, ...) {
	len <- length(unique(difftime(end(x), start(x))))
	return(len == 1)
} )
# tous les intervalles sont de même durée et espacés d'une même période
setMethod (f='regular', signature='TimeIntervalDataFrame',
definition=function(x, ...) {
	len <- length(unique(difftime(start(x)[-1], start(x)[-nrow(x)])))
	return(homogeneous(x) & length(len) == 1)
} )
# entre le début du premier intervalle et la fin du dernier, il
# n'y a pas de 'trous' ET il n'y a pas de superposition entre deux
# intervalles
setMethod (f='continuous', signature='TimeIntervalDataFrame',
definition=function(x, ...) {
	start <- start(x)
	end <- end(x)
	ordre <- order (start, end)
	start <- start[ordre]
	end <- end[ordre]
	return(all (start[-1] == end[-nrow(x)]) )
} )

setMethod (f='continuous<-', signature='TimeIntervalDataFrame',
definition=function(x, value) {
	if (!value) return (x)
	if (overlapping (x) ) stop ("Can't make continuous a 'TimeIntervalDataFrame' with overlapping intervals.")
	if (continuous (x) ) return (x)

	data <-  as.data.frame (matrix(NA, nrow=nrow(x)-1, ncol=ncol(x) ) )
	names (data) <- names (x)

	complementaire <- new ('TimeIntervalDataFrame',
			       start=end(x)[-nrow(x)], end=start(x)[-1],
			       timezone=x@timezone,
			       data= data)

	return (merge (x, complementaire, tz=timezone(x) ) )
} )

# certains intervalles se superposent-ils ?
setMethod (f='overlapping', signature='TimeIntervalDataFrame',
definition=function(x, ...) {
	ol <- .C ('overlapping_timeintervaldf',
		  as.integer(start(x)),
		  as.integer(end(x)),
		  as.integer(nrow(x)),
		  ol=integer(1),
		  NAOK=FALSE, PACKAGE='timetools')$ol
	return (ol == 1)
} )

# transformateur de classe
#-------------------------
setAs ('TimeIntervalDataFrame', 'data.frame',
       function(from) data.frame (start=start(from), end=end(from), from@data) )
as.data.frame.TimeIntervalDataFrame <- function (x, row.names=NULL, optional=FALSE, include.dates=FALSE, ...) {
	if (include.dates)
		return (data.frame (start=start(x), end=end(x), x@data) ) else
		return (x@data)
}

as.TimeInstantDataFrame.TimeIntervalDataFrame <- function(from, cursor=NULL, ...) {
	if (is.null (cursor) ) cursor <- 0.5
	if (cursor > 1 || cursor < 0) warning ("For a standard use, cursor should be between 0 and 1.")
	if(nrow(from) == 0)
		instant <- character(0) else 
		instant <- mapply (function(x, y, wx, wy)
				weighted.mean (c(x, y), c(wx, wy), na.rm=TRUE),
			   start(from), end(from), 1-cursor, cursor)
	instant <- as.POSIXct(instant, origin=timetools::origin)
	to <- new ('TimeInstantDataFrame', instant=instant,
		   timezone=timezone(from), data=from@data)
	validObject(to)
	return (to)
}

as.SubtimeDataFrame.TimeIntervalDataFrame <-
	function(x, unit, of, FUN=NULL, cursor=NULL, ...)
	as.SubtimeDataFrame (as.TimeInstantDataFrame (x, cursor), unit, of, FUN, ...)

