timeDefaults <- function()
{
  x <- timeDateOptions()[ c(
  "time.month.name", "time.month.abb",  "time.day.name",
   "time.day.abb", "time.am.pm", "time.century", "time.zone" ) ]
  list( month.name = x$time.month.name, month.abb = x$time.month.abb,
        day.name = x$time.day.name, day.abb = x$time.day.abb,
        am.pm = x$time.am.pm, century = x$time.century, zone = x$time.zone )
}

timeDateFormatChoose <- function( ms, zone )
{
  ## undocumented function to choose a time/date format
  defform <- timeDateOptions("time.out.format")[[1]]
  notimeform <- timeDateOptions("time.out.format.notime")[[1]]
  ## see if there is a need for time of day
  lst <- timeZoneList()
  ms <- ms[!is.na(ms)]
  hasnotime <- length(ms) && all(ms == 0) &&
      identical(lst[[zone]], lst[["GMT"]])
  if(!hasnotime || !length( notimeform ))
    return( as(defform, "character"))
  return( as(notimeform, "character") )
}

timeDate <- function( charvec, in.format, format, zone, julian, ms,
		  in.origin = c( month = 1, day = 1, year = 1960 ))
{
# Creation function for timeDate objects
# can pass in zero arguments, but if you pass in any,
# either charvec (vector of character formatted times), julian (julian
# days),  or ms (milliseconds) must be supplied
# in.origin can be used with the julian argument to give the date origin
# of the passed in julian days -- integer vector c( month, day, year )

  if( missing( charvec ) && missing( in.format ) &&
     missing( format ) && missing( zone ) && missing( julian )
     && missing( ms ))
  {
    out <- new( "timeDate" )
    return(out)
  }

  # try to use first argument as charvec if it's there
  if( missing( zone )) zone <- timeDateOptions( "time.zone" )[[1]]

  if( !missing( charvec ))
  {
    if( !missing( julian ) || !missing( ms ))
      warning( paste( "Ignoring julian, and/or ms arguments to function",
		      "time, since charvec argument was given" ))

    # test for character vs old-style ts call
    if( !is( charvec, "character" ))
      stop( paste( "First (charvec) argument to function time must",
		   "be of type character" ))

    # get missing arguments
    if( missing( in.format ))
      in.format <- timeDateOptions( "time.in.format" )[[1]]

    # read from the character strings
    defaults <- timeDefaults()
    defaults$zone <- zone

    obj <- .Call( "time_from_string",
                 as(charvec, "character"),
                 as(in.format, "character"),
                 defaults, timeZoneList())
    if( class( obj ) != "timeDate" )
      stop( "Unknown error in calling C function time_from_string" )

  } else {

    # check and take care of missing ms and julian args
    if(missing( ms ))
    {
      if( missing( julian ))
      stop( "time function requires either julian or ms argument" )
      ms <- rep( 0, length( julian ))
    }
    if( missing( julian ))
      julian <- rep( 0, length( ms ))

    # add in the origin
    origin <- .Call( "time_from_month_day_year",
                    as.integer(in.origin[["month"]]),
                    as.integer(in.origin[["day"]]),
                    as.integer(in.origin[["year"]]))
    julian <- julian + groupVecColumn( origin, "julian.day" )

    if( length( ms ) != length( julian ))
      stop( "julian and ms arguments to time must be same length" )

    # create new time object
    # treat julian as numeric if ms is all 0

    if( all( ms == 0, na.rm = TRUE )) {
      if(any(is.na(ms)))
        julian <- julian + ms # copy NA's from ms to julian (non-NA's are 0's)
      obj <- as( julian, "timeDate" )
    }
    else {
      obj <- new( "timeDate" )

      # put in correct data
      # take into account possibility of shorter args
      lj <- length( julian )
      lm <- length( ms )
      if( lj > lm )
      {
	if(( lj %% lm ) != 0 )
	  stop( "Julian and ms arguments have incompatible lengths" )
	ms <- rep( ms, length = lj )
      } else if( lm > lj )
      {
	if(( lm %% lj ) != 0 )
	  stop( "Julian and ms arguments have incompatible lengths" )
	julian <- rep( julian, length = lm )
      }

      obj@columns <- list( as( julian, "integer" ), as( ms, "integer" ))
    }
  }

  # put on rest of input args
  obj@time.zone <- as( zone, "character" )
  if( missing( format ))
    format <- timeDateFormatChoose(obj@columns[[2]], obj@time.zone)
  obj@format <- as( format, "character" )

  obj
}

timeCalendar <- function( m=NULL, d=NULL, y=NULL, h=NULL, min=NULL,
			   s=NULL, ms=NULL, format, zone )
{
# function to create a timeDate object from date as month, day, year and
# time of day as hours, minutes, seconds, milliseconds
#
# the format and zone arguments are optional; other arguments, if
# present, must be of compatible lengths

  lens <- c( length( m ), length( d ), length( y ), length( h ),
	     length( min ), length( s ), length( ms ))
  datalen <- max( lens )
  if( datalen < 1 )
  {
    # no args
    ret <- timeDate()
    if( !missing( format ))
      ret@format <- format
    if( !missing( zone ))
      ret@time.zone <- zone
    return( ret )
  }

  if( any(( datalen %% lens[lens > 0] ) != 0 ))
    stop( "Arguments have incompatible lengths -- must be multiples" )

  # make all arguments the same length
  if( lens[1] == 0 )
    m <- 1
  if( lens[2] == 0 )
    d <- 1
  if( lens[3] == 0 )
    y <- 1960
  if( lens[4] == 0 )
    h <- 0
  if( lens[5] == 0 )
    min <- 0
  if( lens[6] == 0 )
    s <- 0
  if( lens[7] == 0 )
    ms <- 0

  m <- rep( m, length = datalen )
  d <- rep( d, length = datalen )
  y <- rep( y, length = datalen )
  h <- rep( h, length = datalen )
  min <- rep( min, length = datalen )
  s <- rep( s, length = datalen )
  ms <- rep( ms, length = datalen )

  # create the julian days part of the object

  daytimes <- .Call( "time_from_month_day_year",
                    as.integer(m),  as.integer(d),  as.integer(y))

  # create the ms part of the object
  mstimes <- .Call( "time_from_hour_min_sec",  as.integer(h),  as.integer(min),
                    as.integer(s),  as.integer(ms))

  # merge the two times vectors

  groupVecColumn( daytimes, "milliseconds" ) <-
    groupVecColumn( mstimes, "milliseconds" )

  # convert time zones

  if( missing( zone)) zone <- timeDateOptions( "time.zone" )[[1]]

  daytimes <- .Call( "time_to_zone", daytimes, zone, timeZoneList())

  # put the format on

  daytimes@time.zone <- as( zone, "character" )
  if( missing( format ))
    format <- timeDateFormatChoose(as.integer(daytimes@columns[[2]]),
                                   daytimes@time.zone)
  daytimes@format <- as( format, "character" )

  daytimes
}

"timeZoneConvert" <- 
function(x, zone)
{
	# function to convert a time to a local zone
	# NOTE: NORMALLY YOU SHOULD JUST SET THE time.zone SLOT TO THE
	# DESIRED ZONE!!  This function is only useful when the time
	# object was created from data that was from your local time
	# zone but the creation functions weren't told it wasn't GMT.
	#
	# Internally, all times are stored in GMT, and so once a time
	# has been correctly created, zone conversion means just setting
	# the time zone slot so output will look correct.
	if(x@time.zone == "GMT") {
		# originally this function was only intended to work when
		# x was made in GMT time zone, since that was the default one.
		# We continue to use fast algorithm in that case.
		result <- .Call("time_to_zone", x, zone, timeZoneList())
		result@time.zone <- zone
		result@format <- x@format
		result
	}
	else {
		format.orig <- x@format
		x@format <- "%02m/%02d/%Y %02H:%02M:%02S.%03N"
		strings <- paste(as(x, "character"), zone)
		result <- timeDate(strings, format = format.orig, in.format = 
			"%m/%d/%y %H:%M:%S.%N %Z")
		result@time.zone <- zone
		result
	}
}

setMethod("timeConvert", "timeDate", 
    function( x, to.zone, from.zone ) {
	# Converts to x to a new time zone
	if( !missing( from.zone ))
		warning( "For timeDate objects, from.zone is taken from the time zone stored on the object" )
	x@time.zone <- to.zone
	x
} )


setAs( "timeDate", "character",
      function( from ) .Call( "time_to_string",
                             from,
                             timeDefaults(), timeZoneList())
      )

setAs( "character", "timeDate",
  function( from ) timeDate( from ))

setAs( "timeDate", "numeric",
      function( from ) .Call( "time_to_numeric", from)
      )

setAs( "timeDate", "integer", function( from )
      groupVecColumn( from, "julian.day" ))

setAs( "numeric", "timeDate",
      function( from )
      {
	out <- .Call( "time_from_numeric", as.double(from), "timeDate")
	out@time.zone <- as( timeDateOptions( "time.zone" )[[1]], "character" )
	out@format <- timeDateFormatChoose(out@columns[[2]], out@time.zone)
	out
      }
      )

setAs( "Date", "timeDate",
      function( from )
      {
	from <- as.numeric( from )
        val <- timeDate( julian = from )
        return(val)
      })

setAs("POSIXlt", "timeDate", function(from) {
    val <- timeCalendar(
        m=from$mon+1,
        d=from$mday,
        y=from$year + 1900,
        h=from$hour,
        min=from$min,
        s=from$sec
        )
    tzone <- attr(from, "tzone")
    if(!is.null(tzone)) {
        if(length(tzone) == 1) {
            tz <- tzone[1]
        } else {
            tz <- paste(tzone[2], tzone[3], sep="/")
        }
        val <- timeZoneConvert(val, tz)
    }
    val
    }
)

setAs("POSIXct", "timeDate", function(from) {
    as(as.POSIXlt(from), "timeDate")
    }
)

setMethod( "format", signature(x="timeDate"),
          function(x, ...)
          {
            as(x, "character")
          })
            
       
setMethod( "show", "timeDate", function( object )
{
# show method for time objects
  if( length(object ))
  {
    print.default( as( object, "character" ), quote=FALSE)
  }
  else
    cat( "timeDate()\n")
})

setMethod( "summary", "timeDate", function( object, ... )
{
  nas <- is.na( object )
  tmp <- as( object[!nas], "numeric")
  ret <- as( quantile( tmp, c( 0, .25, .5, .75, 1 )), "timeDate" )
  ret <- c( ret, as( tmp, "timeDate" ))
  ret@format = object@format
  ret@time.zone = object@time.zone
  ret = as( ret, "character")
  ret <- ret[c(1,2,3,6,4,5)]
  names( ret ) <- c( "Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max" )
  if( any( nas ))
    ret[ "NAs" ] <- sum( nas )

  ret <- matrix( ret, nrow = 1, dimnames = list( "", names( ret )))
  oldClass( ret ) <- "table"
  ret
})


setMethod( "[", signature( x = "timeDate", i = "ANY" ),
  function(x, i, ..., drop = TRUE )
  {

    ## handle timeEvent subscripting
    if( is(i, "timeEvent")) {
	nevent <- length(i)
	idx <- logical(length(x))
	for (k in 1:nevent){
		idx <- idx | ((i@columns[[1]][k] <= x) & (i@columns[[2]][k] >= x))
	}
	i <- idx
     }
    x@columns <- lapply( x@columns, "[", i, drop=FALSE )
    x
  })


setReplaceMethod( "[",
  signature( x = "timeDate",  i="ANY", j="ANY", value = "timeDate" ),
  function( x, i, j, ..., value )
  {

    ## handle timeEvent subscripting
    if( is(i, "timeEvent")) {
	nevent <- length(i)
	idx <- logical(length(x))
	for (k in 1:nevent){
		idx <- idx | ((i@columns[[1]][k] <= x) & (i@columns[[2]][k] >= x))
	}
	i <- idx
      }
    x[i] <- value@columns
    x
  })

setMethod( "sort.list", signature( x = "positionsCalendar" ),
function( x, partial = NULL,  na.last = TRUE, decreasing = FALSE,
         method = c("shell", "quick", "radix"))
  sort.list(as( x, "numeric" ), partial, na.last, decreasing, match.arg(method)))

setMethod( "sort", signature( x = "positionsCalendar" ),  
           function( x, decreasing = FALSE, ...) 
         {
	     sl <- sort.list( x, decreasing = decreasing, ... )
	     x[ sl ]
	   })

setMethod("mdy", signature( x = "positionsCalendar" ),
   function(x)
   {
     # return a length 3 list with month, day, year of each element
     obj <- .Call( "time_to_month_day_year", x,
		 timeZoneList())
     if( length( obj ) != 3 )
       stop( "Unknown problem in C function time_to_month_day_year" )
     data.frame( month = obj[[1]], day = obj[[2]], year = obj[[3]] )
   }
)

setMethod("hms", signature( x = "positionsCalendar" ),
   function(x)
   {
     # return a length 4 list with hour, minute, second, ms of each value
     obj <- .Call( "time_to_hour_min_sec", x,
		 timeZoneList())
     if( length( obj ) != 4 )
       stop( "Unknown problem in C function time_to_hour_min_sec" )
     data.frame( hour = obj[[1]], minute = obj[[2]], second = obj[[3]],
		 ms = obj[[4]])
   }
)

setMethod( "wdydy", "positionsCalendar",
  function( x )
  {
     wd <- .Call( "time_to_weekday", x,
		 timeZoneList())
     d <- .Call( "time_to_year_day", x,
		 timeZoneList())
     if( length( d ) != 2 )
      stop( "Unknown problem in C function time_to_year_day" )
     data.frame( weekday = wd, yearday = d[[2]], year = d[[1]] )
  }
)

setMethod( "days", "positionsCalendar",
   function( x )
   {
     # return the day of the month of each element as a factor
     ordered( paste( mdy( x )$day ), paste( 1:31 ))
   })

setMethod( "weekdays", "positionsCalendar",
   function( x, abbreviate=TRUE )
   {
     # return the weekday of each element
     d <- wdydy( x )$weekday + 1
     lbl <- if( abbreviate ) timeDateOptions( "time.day.abb" )[[1]]
            else timeDateOptions( "time.day.name" )[[1]]
     lbl <- as( lbl, "character" )
     d <- lbl[ d ]
     ordered( d, levels = lbl, labels = lbl )
   })

setMethod( "months", "positionsCalendar",
   function( x, abbreviate = TRUE )
   {
     # return the month of each element
     d <- mdy(x)$month
     lbl <- if( abbreviate ) timeDateOptions( "time.month.abb" )[[1]]
            else timeDateOptions( "time.month.name" )[[1]]
     lbl <- as( lbl, "character" )
     d <- lbl[ d ]
     ordered( d, levels = lbl, labels = lbl )
   })

setMethod( "quarters", "positionsCalendar",
   function( x, abbreviate = TRUE )
   {
     # return the quarter of each element
     d <- mdy(x)$month
     lbl <- if(abbreviate) c("1Q", "2Q", "3Q", "4Q") else c("I", "II", "III", "IV")
     d <- lbl[ (d-1) %/% 3 + 1 ]
     ordered( d, levels = lbl, labels = lbl )
   })

setMethod( "years", "positionsCalendar",
   function( x )
   {
     # return the year of each element as an ordered factor
     ordered( mdy(x)$year )
   })


setMethod("hours", signature( x = "positionsCalendar" ),
   function(x)
  {
    # return the hour of each element
    hms(x)$hour
  }
)

setMethod("minutes", signature( x = "positionsCalendar" ),
   function(x)
  {
    # return the minute of each element
    hms(x)$minute
  }
)

setMethod("seconds", signature( x = "positionsCalendar" ),
   function(x)
  {
    # return the seconds of each element as numeric, including
    # fractional seconds
    tms <- hms(x)
    tms$second + ( tms$ms / 1000.0 )
  }
)

setMethod("yeardays", signature( x = "positionsCalendar" ),
   function(x)
{
  # return the day of the year, 1 - 366, of each element
  wdydy(x)$yearday
}
)

setMethod( "Summary", signature( x = "positionsCalendar" ),
	  function(x, ..., na.rm = FALSE)
	  callGeneric( as.numeric(x), ..., na.rm = na.rm ))

setMethod( "Compare", signature( e1 = "positionsCalendar", e2 = "positionsCalendar" ),
	  function( e1, e2 )
	    callGeneric( as( e1, "numeric" ), as( e2, "numeric" )),
	  )

setMethod( "Math", "positionsCalendar",
	  function( x ) callGeneric(as.numeric(x)))

setMethod( "timeFloor", "positionsCalendar",
	  function( x )
	  {
	    # floor method for time objects
	    ret <- .Call( "time_floor", x, timeZoneList() )
	    ret@format <- x@format
	    ret@time.zone <- x@time.zone
	    ret
	  })

setMethod( "timeCeiling", "positionsCalendar",
	  function( x )
	  {
	    # ceiling method for time objects
	    ret <- .Call( "time_ceiling", x, timeZoneList() )
	    ret@format <- x@format
	    ret@time.zone <- x@time.zone
	    ret
	  })

setMethod( "timeTrunc", "positionsCalendar",
	   function( x ) timeFloor( x ))

setMethod( "floor", "positionsCalendar",
	  function( x )
	  {
	    # floor method for time objects
	    ret <- .Call( "time_floor", x, timeZoneList())
	    ret@format <- x@format
	    ret@time.zone <- x@time.zone
	    ret
	  })

setMethod( "ceiling", "positionsCalendar",
	  function( x )
	  {
	    # ceiling method for time objects
	    ret <- .Call( "time_ceiling", x, timeZoneList())
	    ret@format <- x@format
	    ret@time.zone <- x@time.zone
	    ret
	  })

setMethod( "trunc", "positionsCalendar",
	   function( x ) floor( x ))


setMethod( "Math2", signature( x = "positionsCalendar"),
	  function( x, digits ) callGeneric( as.numeric(x), digits ))

setMethod( "Summary", signature( x = "positionsCalendar" ),
	  function(x, ..., na.rm = FALSE)
	  callGeneric( as(x, "numeric"), ..., na.rm = na.rm ))

setMethod( "Ops", signature( e1 = "positionsCalendar", e2 = "ANY" ),
	  function( e1, e2 = NULL )
	  callGeneric( as(e1, "numeric"), e2))

setMethod( "Ops", signature( e1 = "ANY", e2 = "positionsCalendar" ),
	  function( e1, e2 = NULL )
	  callGeneric( e1, as(e2, "numeric")))

setMethod( "Ops", signature( e1 = "positionsCalendar", e2 = "positionsCalendar" ),
	  function( e1, e2 = NULL )
	  callGeneric( as(e1, "numeric"), as(e2, "numeric")))

##setMethod("factor", signature( x = "positionsCalendar" ), 
##          function(x, levels, labels = levels, exclude = NA, 
##                   ordered = is.ordered(x))
##  {
##    x <- as(x, "character")
##    callGeneric(x, levels, labels, exclude, ordered)
##  })

setMethod("shiftPositions", signature(x="vector"),
   function(x, k = 1)
   {
     ## default method for vectors
     if(length(x) <= 1) stop("Cannot shift length 0 or 1 vector")
     diffs <- as(x[-1] - x[ - length(x)], "numeric")
     tol <- timeDateOptions("sequence.tol")[[1]]
     if(any(abs(diffs - diffs[1]) > tol))
       stop("Cannot shift irregular vector")
     x + k * diffs[1]
})

setMethod("shiftPositions", "timeDate", 
function(x, k = 1)
{
        if(k != round(k)) {
          k <- round(k)
          warning("k is not an integer")
        }
        len <- length(x)
        if(len <= 1)
          stop("Cannot shift length 0 or 1 vector")
        dfs <- diff(x)
        numdfs <- as(dfs, "numeric")
        tol <- timeDateOptions("sequence.tol")[[1]]
        if(all(abs(numdfs - numdfs[1]) <= tol)) {
          diffToUse <- k * dfs[1]
        }
        else {
          ## Maybe it's a monthly or yearly sequence
          ## if hours/minutes/seconds all same
          hms1 <- hms(x)
          if(any(c(hms1$hour[-1] - hms1$hour[ - len], hms1$minute[-1] -
                   hms1$minute[ - len], hms1$second[-1] -
                   hms1$second[ - len], hms1$ms[-1] - hms1$ms[ - len])))
            stop("cannot shift irregular vector")
          ## check if all dates are month end dates or same dates
          mdy1 <- mdy(x)
          if(!all(is.monthend(x)) && any(mdy1$day[-1] - mdy1$day[ - len]) )
            stop("cannot shift irregular vector")
          mdiffs <- mdy1$month[-1] - mdy1$month[ - len] + 12 *
            (mdy1$year[-1] - mdy1$year[ - len])
          if(any(mdiffs - mdiffs[1] != 0))
            stop("cannot shift irregular vector")
                diffToUse <- timeRelative(by = "months", k.by = mdiffs[1] * k)
        }
        x + diffToUse
 })

setMethod("as.character", "timeDate", function(x, ...) as(x, "character"))
