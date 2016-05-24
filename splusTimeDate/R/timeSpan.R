
timeSpan <- function( charvec, in.format, format, julian, ms )
{
# creation function for timeSpan objects
# can pass in zero arguments, but if you pass in any,
# either charvec (vector of character formatted time spans), julian (julian
# days),  or ms (milliseconds) must be supplied
  if( missing( charvec ) && missing( in.format ) &&
     missing( format ) && missing( julian )
     && missing( ms ))
  {
    out <- new( "timeSpan" )
    out@format <- as( timeDateOptions("tspan.out.format")[[1]], "character" )
    return(out)
  }


  # try to use first argument as charvec if it's there
  if( !missing( charvec ))
  {
    if( !missing( julian ) || !missing( ms ))
      warning( paste( "Ignoring julian, and/or ms arguments to function",
		      "timeSpan, since charvec argument was given" ))

    # get missing arguments
    if( missing( in.format ))
      in.format <- timeDateOptions( "tspan.in.format" )[[1]]

    # read from the character strings
    obj <- .Call( "tspan_from_string", charvec, in.format)
    if( class( obj ) != "timeSpan" )
      stop( "Unknown error in calling C function tspan_from_string" )

  } else {

    # check and take care of missing ms and julian args

    if(missing( ms ))
    {
      if( missing( julian ))
      stop( "timeSpan function requires either julian or ms argument" )
      ms <- rep( 0, length( julian ))
    }
    if( missing( julian ))
      julian <- rep( 0, length( ms ))


    # treat julian as numeric if ms is all 0

    if( all( !is.na(ms) & ms == 0 ))
      obj <- as( julian, "timeSpan" )
    else
    {
      # create new time span object
      obj <- new( "timeSpan" )

      # take into account possibility of shorter args
      lj <- length( julian )
      lm <- length( ms )
      if( lj > lm ) {
        if(( lj %% lm ) != 0 )
          stop( "Julian and ms arguments have incompatible lengths" )
        ms <- rep( ms, length = lj )
      } else if( lm > lj ) {
        if(( lm %% lj ) != 0 )
          stop( "Julian and ms arguments have incompatible lengths" )
        julian <- rep( julian, length = lm )
      }
      # put in correct data
      jul.int <- as(julian, "integer")
      ms.int <- as(ms, "integer")
      ex <- !is.na(ms.int) & ms.int>=86400000
      if (any(ex)) {
        jul.int[ex] <- jul.int[ex] + (ms.int[ex] %/% 86400000)
        ms.int[ex] <- ms.int[ex] %% 86400000
      }
      ex <- !is.na(ms.int) & ms.int<=-86400000
      if (any(ex)) {
        jul.int[ex] <- jul.int[ex] - ((-ms.int[ex]) %/% 86400000)
        ms.int[ex] <- - ((-ms.int[ex]) %% 86400000)
      }
      obj@columns <- list(jul.int, ms.int)
    }
  }

  # put on format
  if( missing( format )) format <- timeDateOptions("tspan.out.format")[[1]]
  obj@format <- as( format, "character" )

  obj
}

setAs( "timeSpan", "character",
      function( from ) .Call( "tspan_to_string", from)
      )

setAs( "character", "timeSpan",
  function( from ) timeSpan( charvec = from ))

setAs( "timeSpan", "numeric",
      function( from ) .Call( "time_to_numeric", from)
      )

setAs( "timeSpan", "integer", function( from )
      groupVecColumn( from, "julian.day" ))

setAs( "numeric", "timeSpan",
      function( from )
      {
	out <- .Call( "time_from_numeric", as.double(from), "timeSpan")
	out@format <- as( timeDateOptions("tspan.out.format")[[1]], "character" )
	out
      })

# setMethod( "format", signature(x="timeDate"),
setMethod( "format", signature(x="timeSpan"),
          function(x, ...)
          {
            as(x, "character")
          })

setMethod( "show", "timeSpan", function( object )
{
# show method for timeSpan objects
  if( length( object ))
  {
     print.default( as( object, "character" ), quote=FALSE)
  }
  else
    cat( "timeSpan()\n" )
})

setMethod( "summary", "timeSpan", function( object, ... )
{
  nas <- is.na( object )
  ret <- as( quantile( object[!nas], c( 0, .25, .5, .75, 1 )), "character" )
  ret <- c( ret, as( mean( object[!nas] ), "character" ))
  ret <- ret[c(1,2,3,6,4,5)]
  names( ret ) <- c( "Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max" )
  if( any( nas ))
    ret[ "NAs" ] <- sum( nas )

  ret <- matrix( ret, nrow = 1, dimnames = list( "", names( ret )))
  oldClass( ret ) <- "table"
  ret
})

setMethod( "sort.list", signature( x = "timeSpan" ),
function( x, partial = NULL,  na.last = TRUE, decreasing = FALSE,
         method = c("shell", "quick", "radix"))
  sort.list( as( x, "numeric" ), partial, na.last, decreasing, match.arg(method)))

setMethod( "sort", signature( x = "timeSpan" ),
           function( x, decreasing = FALSE, ...)
           { 
	     sl <- sort.list( x, decreasing = decreasing, ... )
	     x[ sl ]
	   })

setMethod( "Math", "timeSpan",
	  function( x ) callGeneric(as(x, "numeric")),
	  )

setMethod( "timeFloor", "timeSpan",
	  function( x )
	  {
	    # floor subtracts one to julian day wherever ms < 0
	    where.to.sub <- ( groupVecColumn( x, "milliseconds" ) < 0 )
	    ret.value <- groupVecColumn( x, "julian.day" )
	    ret.value[ where.to.sub ] <- ret.value[ where.to.sub ] - 1
	    groupVecData( x ) <- list( as.integer(ret.value),
				  rep( 0L, length( ret.value )))
	    x
	  })

setMethod( "timeCeiling", "timeSpan",
	  function( x )
	  {
	    # ceiling adds one to julian day wherever ms > 0
	    where.to.add <- ( groupVecColumn( x, "milliseconds" ) > 0 )
	    ret.value <- groupVecColumn( x, "julian.day" )
	    ret.value[ where.to.add ] <- ret.value[ where.to.add ] + 1
	    groupVecData( x ) <- list( as.integer(ret.value),
				  rep( 0L, length( ret.value )))
	    x
	  })

setMethod( "timeTrunc", "timeSpan",
	  function( x )
	  {
	    # we truncate by taking julian day part of the time span
	    ret.value <- groupVecColumn( x, "julian.day" )
	    groupVecData( x ) <- list( as.integer(ret.value),
				  rep( 0L, length( ret.value )))
	    x
	  })

setMethod( "Math2", signature( x = "timeSpan"),
	  function( x, digits ) callGeneric(as(x, "numeric"), digits),
	  )

setMethod( "Summary", signature( x = "timeSpan" ),
	  function(x, ..., na.rm = FALSE)
	  callGeneric( as(x, "numeric"), ..., na.rm = na.rm ))

setMethod( "Ops", signature( e1 = "timeSpan", e2 = "ANY" ),
	  function( e1, e2 = NULL )
	  callGeneric( as(e1, "numeric"), e2))

setMethod( "Ops", signature( e1="ANY", e2 = "timeSpan" ),
          function( e1, e2 = NULL )
	  callGeneric( e1, as(e2, "numeric")))

setMethod( "Ops", signature( e1="timeSpan", e2 = "timeSpan" ),
          function( e1, e2 = NULL )
	  callGeneric( as(e1, "numeric"), as(e2, "numeric")))

setMethod("hms", signature(x="timeSpan"),
    function(x) {
        class(x) <- "timeDate"
        x@time.zone <- "GMT"
        obj <- .Call("time_to_hour_min_sec", x, timeZoneList())
        if(length(obj) != 4) {
            stop("Unknown problem in C function time_to_hour_min_sec")
        }
        data.frame(hour = obj[[1]], minute = obj[[2]], second = obj[[3]],
            ms = obj[[4]])
     }
)
        
setMethod("hours", signature(x="timeSpan"),
    function(x) {
        class(x) <- "timeDate"
        x@time.zone <- "GMT"
        callGeneric()
    }
)

setMethod("minutes", signature(x="timeSpan"),
    function(x) {
        class(x) <- "timeDate"
        x@time.zone <- "GMT"
        callGeneric()
    }
)

setMethod("seconds", signature(x="timeSpan"),
    function(x) {
        class(x) <- "timeDate"
        x@time.zone <- "GMT"
        callGeneric()
    }
)

setMethod("as.character", "timeSpan", function(x, ...) as(x, "character"))

##setMethod("factor", signature( x = "timeSpan" ), 
##          function(x, levels, labels = levels, exclude = NA, 
##                   ordered = is.ordered(x))
##  {
##    x <- as(x, "character")
##    callGeneric(x, levels, labels, exclude, ordered)
##  })

