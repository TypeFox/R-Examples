timeSequence <- function( from, to, by, length.out, k.by=1, align.by=FALSE,
                         extend=FALSE, week.align=NULL, holidays=timeDate(),
			 exceptions, additions, format, zone )
{
  if( missing( from ) && missing( to ) && missing( by ) &&
      missing( length.out ) && missing( exceptions ) &&
      missing( additions ) && missing( format ) && 
      missing( zone ))
  {
    out <- new( "timeSequence" )
    out@format <- as( timeDateOptions( "time.out.format" )[[1]], "character")
    out@time.zone <- as( timeDateOptions( "time.zone" )[[1]], "character")
    return(out)
  }

  if(missing(zone)) {
    zone <- as(timeDateOptions("time.zone")[[1]], "character")
  }
  else {
    time.zone <- as(timeDateOptions(time.zone = zone)[[1]], "character")
    on.exit(timeDateOptions(time.zone = time.zone))
  }

  if( missing( from )) from <- timeDate() 
    else from <- as( from, "timeDate")

  if( missing( to )) to <- timeDate()
    else to <- as( to, "timeDate")

  if( missing( length.out )) length.out <- numeric(0)
  length.out <- as( length.out, "integer" )

  if( missing( by )) by <- timeSpan()

  if( is( by, "character" ))
  {
    # check sign of k.by against from and to
    if( length( from ) && length( to ) && 
        ( k.by * as( to - from, "numeric" ) < 0 ))
      stop( "Incompatible values of from and to" )

    # align the from and to times
    if( align.by )
    {
      if (by == "quarters") {
        if (4 %% k.by != 0) stop("when using by=\"", by, "\", k.by must be an even divisor of 4")
      } else if (by == "months") {
        if (12 %% k.by != 0) stop("when using by=\"", by, "\", k.by must be an even divisor of 12")
      } else if (by == "days" || by == "bizdays" || by == "weekdays" ) {
        if (k.by <= 0 || k.by > 30) stop("when using by=\"", by, "\", k.by must be in the range: 1:30")
      } else if (by == "hours") {
        if (24 %% k.by != 0) stop("when using by=\"", by, "\", k.by must be an even divisor of 24")
      } else if (by == "minutes" || by == "seconds") {
        if (60 %% k.by != 0) stop("when using by=\"", by, "\", k.by must be an even divisor of 60")
      }
    
      reversed <- length( from ) && length( to ) && ( from > to )
      if( length( from ))
      {
	if( extend == reversed )
	  from <- timeAlign( from, by, k.by, 1, week.align )
	else
	  from <- timeAlign( from, by, k.by, -1, week.align )
      }
      if( length( to ))
      {
	# special case for "weeks" if they didn't specify week.align,
	# take it from the weekday of the from, if it's there...
	if( length( from ) && is.null( week.align ) && 
		by == "weeks" ) {
		tmp.week.align <- wdydy( from )$weekday
	} else {
		tmp.week.align <- week.align
	}
	if( extend == reversed )
	  to <- timeAlign( to, by, k.by, -1, tmp.week.align )
	else
	  to <- timeAlign( to, by, k.by, 1, tmp.week.align )
      }
    }

    # special case for weekdays/business days -- move "from" and "to"
    # to the next week/business day just in case it is not already
    # but don't align to midnight (if align.by = TRUE, we've already done this)
    # also for weeks, with a given week.align, move to correct weekday
    if(!align.by && ( by == "bizdays" || by == "weekdays" || 
	( by == "weeks" && !is.null(week.align) ))) {
       reversed <- length( from ) && length( to ) && ( from > to )
       if( by == "bizdays" ) 
	   str <- "biz"
       else 
	   str <- "wkd"

       if( by == "weeks" ) {
          tonext <- timeRelative( by="weeks", k.by=1, week.day = week.align )
          toprev <- timeRelative( by="weeks", k.by=-1, week.day = week.align )
       } else {
          tonext <- timeRelative( paste("-1", str, " +1", str, sep = ""),
              holidays. = holidays )
          toprev <- timeRelative(paste("+1", str, " -1", str, sep = ""),
              holidays.=holidays )
       }

      if( reversed == extend ) {
	  from <- from + tonext
	  to <- to + toprev
       } else {
	  from <- from + toprev
	  to <- to + tonext
       }
    }

    # now see if alignment has nullified us 
    if( length( from ) && length( to ) && 
        ( k.by * as( to - from, "numeric" ) < 0 ))
    {
      to <- timeDate()
      length.out <- 0L
    }

    # convert to a timeRelative object
    by <- timeRelative( by=by, k.by=k.by, align.by=FALSE,
			 week.day=week.align, holidays.=holidays)
  } else if( !missing( k.by ) || !missing( align.by ) || 
	     !missing( extend ) || !missing( week.align ))
    warning( "k.by, align.by, extend, and week.align are ignored if by is not a character string" )

  if( !is( by, "timeInterval" ))
    by <- as( by, "timeSpan" )

  lens <- c( length(from), length(to), length(by), length(length.out))
  if( any( lens > 1 ) || ( sum(lens) != 3 ))
    stop( paste( "Exactly three of from, to, by, length.out must be specified",
		"and they must each be of length 1" ))
  if( any( c( is.na( from ), is.na( to ), is.na(by), is.na(length.out))))
    stop( "NA values not allowed in sequences" )

  ret <- new( "timeSequence", from = from, to = to, by = by, 
	      length = as.integer(length.out), time.zone = zone )

  if( !missing( exceptions ))
    ret@exceptions <- as( exceptions, "timeEvent" )
  if( !missing( additions ))
    ret@additions <- as( additions, "positionsCalendar" )
  ret@time.zone <- as( zone, "character" )
  if( missing( format )) {
    # to choose format, we need to see if we have any milliseconds
    # on the dateTime objects
    tofromms <- c( from@columns[[2]], to@columns[[2]] )
    if( is( by, "timeSpan" ) ) {
	    tofromms <- c( tofromms, by@columns[[2]] )
    } else {
        if( length(length.out) && !length(by)) {
	       # if they give to/from/length, assume we'll get to fractional days
	        tofromms <- c( tofromms, 1 )
        } else {
            ## we should have by and either from or to
	        if( length( from )) {
                nextone <- from + by
	        } else {
                nextone <- to - by
            }
	        tofromms <- c( tofromms, nextone@columns[[2]] )
        }
    }
    format <- timeDateFormatChoose(tofromms, ret@time.zone)
  }
  ret@format <- as( format, "character" )
  ret
}

setMethod( "format", signature(x="timeDate"),
          function(x, ...)
          {
            as(x, "character")
          })

setMethod( "show", "timeSequence", function( object )
{
  if( !is.null( object@from ) && length(object@from)>0 )
  {
    tmp <- object@from
    tmp@time.zone <- object@time.zone
    tmp@format <- object@format
    tmp <- as( tmp, "character" )
    cat( "from:  ", tmp, "\n" )
  }
  if( !is.null( object@to ) && length(object@to)>0 )
  {
    tmp <- object@to
    tmp@time.zone <- object@time.zone
    tmp@format <- object@format
    tmp <- as( tmp, "character" )
    cat( "to:    ", tmp, "\n" )
  }
  if( !is.null( object@by ) && length(object@by)>0 )
    cat( "by:    ", as( object@by, "character"), "\n" )
  if( !is.null( object@length ) && length(object@length)>0 )
    cat( "length (before exceptions/additions):", object@length, "\n" )
  if( !is.null( object@by )  && length(object@by)>0 &&
     is( object@by, "timeRelative" ) &&
     length( object@by@holidays )>0 )
    cat( "business day holidays:", length(object@by@holidays), "\n" )
  if( length( object@exceptions )>0 )
    cat( "exceptions:", length(object@exceptions), "\n" )
  if( length( object@additions )>0 )
    cat( "additions:", length(object@additions), "\n" )
  dts <- as(as(object, "timeDate"), "character")
  if( length(dts) > 5 )
    dts <- c( dts[1:3], "...", dts[length(dts)])
  invisible(print(dts, quote=FALSE))
})

setMethod( "summary", "timeSequence", 
   function( object, ... )
   {
     vec <- class(object)
     nms <- "Class"

     if( !is.null( object@from ) && length(object@from)>0 )
     {
       tmp <- object@from
       tmp@time.zone <- object@time.zone
       tmp@format <- object@format
       tmp <- as( tmp, "character" )
       vec <- c( vec, tmp )
       nms <- c( nms, "From" )
     }

     if( !is.null( object@to ) && length(object@to)>0 )
     {
       tmp <- object@to
       tmp@time.zone <- object@time.zone
       tmp@format <- object@format
       tmp <- as( tmp, "character" )
       vec <- c( vec, tmp )
       nms <- c( nms, "To" )
     }

     if( !is.null( object@by ) && length(object@by)>0 )
     {
       vec <- c( vec, as( object@by, "character" ))
       nms <- c( nms, "By" )
     }

     if( !is.null( object@length ) && length(object@length)>0 )
     {
       vec <- c( vec, object@length )
       nms <- c( nms, "Length" )
     }

     if( !is.null( object@by ) && length(object@by)>0 &&
        is( object@by, "timeRelative" ) &&
	length( object@by@holidays ))
     {
       vec <- c( vec, length( object@by@holidays ))
       nms <- c( nms, "Biz.Holidays" )
     }

     if( length(object@exceptions))
     {
       vec <- c( vec, length( object@exceptions ))
       nms <- c( nms, "Exceptions" )
     }

     if( length(object@additions))
     {
       vec <- c( vec, length( object@additions ))
       nms <- c( nms, "Additions" )
     }

     ret <- matrix( vec, nrow=1, dimnames=list("", nms ))
     oldClass(ret) <- c( "table", "matrix" )
     ret
   })

setAs( "timeSequence", "timeDate", 
function(from)
  {
    if( is.null( from@by ) || length(from@by)==0 )
      {
        ## we have from, to, and length
        ret <- as( seq( from = as( from@from, "numeric" ),
                       to = as( from@to, "numeric" ),
                       length = from@length ), "timeDate" )
      } else if( is( from@by, "timeSpan" ))
	{
	  by <- as( from@by, "numeric" )
	  if( is.null( from@from) || length(from@from)==0 )
	    ret <- as( seq( to = as( from@to, "numeric" ),
                           by=by, length = from@length ), "timeDate")
	  else if( is.null( from@to ) || length(from@to)==0 )
	    ret <- as( seq( from = as( from@from, "numeric" ),
                           by=by, length = from@length ), "timeDate")
	  else
	    ret <- as( seq( from = as( from@from, "numeric" ),
                           to = as( from@to, "numeric" ),
                           by=by ), "timeDate")
	} else 
    {
      ## if we get here, by is a timeRelative from
      if( length( from@from )==0)
        ret <- rev( .Call( "time_rel_seq",
                          as(from@to, "timeDate"),
                          timeDate(), 
                          as(from@length, "integer"),
                          TRUE,
                          as((-from@by)@Data, "character"),
                          as(sort(from@by@holidays), "timeDate"),
                          timeZoneList()))
      else if( length( from@to )==0)
        ret <- .Call( "time_rel_seq",
                     as(from@from, "timeDate"),
                     timeDate(), 
                     as(from@length, "integer"),
                     TRUE,
                     as(from@by@Data, "character"),
                     as(sort(from@by@holidays), "timeDate"),
                     timeZoneList())
      else
        ret <- .Call( "time_rel_seq",
                     as(from@from, "timeDate"),
                     as(from@to,"timeDate"),
                     0L,
                     FALSE,
                     as(from@by@Data, "character"),
                     as(sort(from@by@holidays), "timeDate"),
                     timeZoneList())
    }
    ## do exceptions and additions
    rever <- (( length( ret ) > 1 ) && ( ret[1] > ret[2] ))
    
    exc <- groupVecColumn( from@exceptions, c( "start", "end" ))

    ## exc is now a list whose first element is the vector of starts
    ## and second is the vector of ends.
    ## so now we have to cut out the times that are in these
    ## intervals
    
    if( length( exc[[1]] ) && length( ret ))
      {
        rng <- range( ret )
        ## assume starts/ends are in order
        starts <- exc[[1]]
        ends <- exc[[2]]
        strt <- min( rng[1], starts[1]) - 1
        endt <- max( rng[2], ends[length(ends)]) + 1
        ## separate ret into bins by the "starts" and "ends"...
        catstrt <- as.numeric(cut( ret, c( strt, exc[[1]], endt ),
                                  right=FALSE)) - 1
        catend <- as.numeric(cut( ret, c( strt, exc[[2]], endt ),
                                 right=TRUE))
        ## ... in order to figure out which are between a pair of start/end
        inperiod <- catstrt == catend
        if( any( inperiod ))
          ret <- ret[ !inperiod ]
      }
    ret <- c( ret, from@additions )
    ret <- unique( sort( ret ))
    if( rever ) ret <- rev( ret )
    ## put on format and zone
    ret@format <- from@format
    ret@time.zone <- from@time.zone
    ret
  }
)
       


setAs( "timeSequence", "character", 
       function(from) as( as( from, "timeDate" ), "character" )
       )

setAs( "timeSequence", "numeric", 
      function( from ) as( as( from, "timeDate" ), "numeric" )
      )

setAs( "timeSequence", "integer", 
      function( from ) as( as( from, "timeDate" ), "integer" )
      )

## modified by Quan Wen on 9/30/04 to detect all regular time series that fall on month ends
## but still there will be problems with dates generated by timeCalendar(), timeSequence() or timeSeq() functions
##	that should fall on month end sequence of dates but not.  

setAs( "timeDate", "timeSequence",
      function(from)
      {
	# convert timeDate from to a sequence, if possible 
	len <- length( from )
	if( !len )
	  return( timeSequence() )
	if( len <= 2 )
	  return( timeSequence( from=from[1], to=from[len], length.out = len ))

	# Try as a regularly-spaced sequence
	diffs <- from[ -1 ] - from[ - len ]
	tol <- timeDateOptions( "sequence.tol" )[[1]]
	if( !any(abs( as(diffs - diffs[1], "numeric" )) > tol ))
	  return( timeSequence( from = from[1], by = diffs[1], 
			 length.out = len ))

	# Maybe it's a monthly or yearly sequence.  First check hours, min
	# sec, etc and days match
	hms1 <- hms( from )
	mdy1 <- mdy( from )
	if( any(c( hms1$hour[-1] - hms1$hour[ -len ],
		     hms1$minute[-1] - hms1$minute[ -len ],
		     hms1$second[-1] - hms1$second[ -len ],
		     hms1$ms[-1] - hms1$ms[ -len ] )))
	  stop( "Cannot detect regular pattern in sequence" )
## check if all dates are month end dates
	if (!all(is.monthend(from))) {
		if( any(mdy1$day[-1] - mdy1$day[ -len ]))
			  stop( "Cannot detect regular pattern in sequence" )			
	}
	mdiffs <- mdy1$month[-1] - mdy1$month[-len] + 
	  12 * ( mdy1$year[-1] - mdy1$year[-len] )
	if( any( mdiffs - mdiffs[1] != 0 ))
	  stop( "Cannot detect regular pattern in sequence" )

	return( timeSequence( from = from[1], length.out = len, 
			      by = "months", k.by = mdiffs[1],
			      format = from@format, 
			      zone=from@time.zone))

      })

setMethod("as.character", "timeSequence", function(x, ...) as(x, "character"))

setMethod( "c", "timeSequence", 
  function(x,...) c(as(x, "timeDate"), ... ))

setMethod( "is.na", "timeSequence",
  function(x) is.na(as(x, "timeDate")))

setMethod( "duplicated", "timeSequence",
  function(x, incomparables=FALSE) 
    duplicated(as(x, "timeDate"), incomparables))

setMethod( "length", "timeSequence",
  function(x) length(as(x, "timeDate")))

setMethod( "[", "timeSequence", 
   function(x, ..., drop = TRUE )
   {
     x <- as( x, "timeDate" )
     callGeneric()
   })

setMethod( "[[", "timeSequence", 
   function(x, ... )
   {
     x <- as( x, "timeDate" )
     callGeneric()
   })

setReplaceMethod( "[[", "timeSequence",
   function(x, ..., value )
   { 
     x <- as( x, "timeDate" )
     x[[...]] <- value
     x
   })

setReplaceMethod( "[", "timeSequence",
   function(x, ..., value )
   { 
     x <- as( x, "timeDate" )
     x[...] <- value
     x
   })


setMethod( "length<-", "timeSequence",
         function( x, value )
       {
	x <- as( x, "timeDate" )
	callGeneric()
       })

"timeSeq" <- 
function(from, to, by = "days", length.out, k.by = 1, align.by = FALSE,
extend = FALSE,	week.align = NULL, holidays, exceptions, additions,
format, zone)
{
  if(missing(from))
    seq.call <- list(to = to, by = by, length.out = length.out)
  else if(missing(to))
    seq.call <- list(from = from, by = by, length.out = length.out)
  else {
    if(!missing(length.out) && !missing(by))
      warning("sequence is over-specified -- ignoring length"
              )
    if(!missing(length.out) && missing(by))
      seq.call <- list(from = from, to = to, length.out = 
                       length.out)
    else seq.call <- list(from = from, to = to, by = by)
  }
  if(!missing(k.by))
    seq.call[["k.by"]] <- k.by
  if(!missing(align.by))
    seq.call[["align.by"]] <- align.by
  if(!missing(extend))
    seq.call[["extend"]] <- extend
  if(!missing(week.align))
    seq.call[["week.align"]] <- week.align
  if(!missing(holidays))
    seq.call[["holidays"]] <- holidays
  if(!missing(exceptions))
    seq.call[["exceptions"]] <- exceptions
  if(!missing(additions))
    seq.call[["additions"]] <- additions
  if(!missing(format))
    seq.call[["format"]] <- format
  else if(!missing(from) && is(from, "timeDate"))
    seq.call[["format"]] <- from@format
  else if(!missing(to) && is(to, "timeDate"))
    seq.call[["format"]] <- to@format
  if(!missing(zone))
    seq.call[["zone"]] <- zone
  else if(!missing(from) && is(from, "timeDate"))
    seq.call[["zone"]] <- from@time.zone
  else if(!missing(to) && is(to, "timeDate"))
    seq.call[["zone"]] <- to@time.zone
  ret <- do.call("timeSequence", seq.call)
  as(ret, "timeDate")
}

setMethod("shiftPositions", signature(x="timeSequence"),
	  function(x, k=1)
	  {
	    if(k != round(k)) {
            k <- round(k)
	        warning("k is not an integer")
	    }
	    xlength <- length(x)
	    xby <- x@by
	    if(!length(xby)) {
		    xby <- (x@to - x@from) / (xlength - 1)
	    }
	    if(k < 0) {
	        x <- timeSequence(to=x[xlength], by=xby,
	            length.out=xlength + abs(k))
		    timeSequence(from=x[1], by=xby, length.out=xlength)
        } else {
            if(k >= xlength) {
                x <- timeSequence(from=x[1], by=xby, length.out=xlength+k)
            }
            timeSequence(from=x[1 + k], by=xby, length.out=xlength)
	    }
	  })
