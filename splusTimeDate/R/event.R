"timeEvent" <- 
function(start., end., IDs)
{
  ret <- new("timeEvent")
  if(missing(start.) && missing(end.) && missing(IDs))
    return(ret)
  if(!is(start., "positionsCalendar"))
    start. <- as(start., "timeDate")
  ## if end is missing, make it 1 day - 1 ms past start
  if(missing(end.)) end. <- start. + (1 - 1/(3600 * 24 * 1000))
  if(!is(end., "positionsCalendar"))
    end. <- as(end., "timeDate")
  if(missing(IDs))
    IDs <- rep("", length(start.))
  if((length(end.) != length(start.)) || (length(end.) != length(IDs)))
    stop("Data lengths must match")
  ret@columns <- list(start., end., IDs)
  ret
}

setMethod( "show", "timeEvent", 
	   function( object )
	   {
	     if( !length( object ))
	     {
	       cat( "timeEvent()\n" )
	       return( NULL )
	     }

	     nms <- as( groupVecColumn( object, "IDs" ), "character" )
	     starts <- as( groupVecColumn( object, "start" ), "character" )
	     ends <- as( groupVecColumn( object, "end" ), "character" )

	     outmat <- cbind( nms, starts, ends )
	     dimnames( outmat ) <- 
	       list( rep( "", length( nms )), c( "ID", "start", "end" ))
	     oldClass( outmat ) <- "table"
	     show( outmat )

	     NULL
	   })

setMethod( "summary", "timeEvent", function( object, ... )
{
  starts <- groupVecColumn( object, "start" )
  nas <- is.na( starts )
  starts <- c( as( range( starts[ !nas ] ), "character"), sum( nas ) )

  ends <- groupVecColumn( object, "end" )
  nas <- is.na( ends )
  ends <- c( as( range( ends[ !nas ] ),  "character"), sum( nas ))

  ret <- cbind( starts, ends )
  dimnames( ret ) <- list( c( "Min", "Max", "NAs" ), c( "Start", "End" ))
		                  
  oldClass( ret ) <- "table" 
  ret
})

setAs( "positionsCalendar", "timeEvent",
       function( from ) timeEvent( from ))

setMethod( "Math", "timeEvent", 
	  function( x ) stop( "Function does not make sense for events" )
          )

setMethod( "Math2", signature( x = "timeEvent"),  
	  function( x, digits ) 
	    stop( "Function does not make sense for events" ))

setMethod( "Summary", signature( x = "timeEvent" ),
	  function(x, ..., na.rm = FALSE) 
	    stop( "Function does not make sense for events" ))

setMethod( "Ops", signature( e1 = "timeEvent" ),
	  function( e1, e2 = NULL )
	    stop( "Function does not make sense for events" ))

setMethod( "Ops", signature( e2 = "timeEvent" ),
	  function( e1, e2 = NULL )
	    stop( "Function does not make sense for events" ))

setMethod( "==", signature( e1 = "timeEvent", e2 = "timeEvent" ),
	  function( e1, e2 = NULL )
{
  dat1 <- e1@columns
  dat2 <- e2@columns
  (dat1[[1]] == dat2[[1]] ) & (dat1[[2]] == dat2[[2]] ) & 
    (dat1[[3]] == dat2[[3]] )
})

setMethod( "initialize", "timeEvent",
          function( .Object, names, classes, columns, ...){
            if(missing(names))
              names <- c( "start", "end", "IDs" )
            if(missing(classes))
              classes <- c( "positionsCalendar", "positionsCalendar", "ANY" )
            if(missing(columns))
              columns <- list( new("timeDate"), new("timeDate"), character() )
            callNextMethod(.Object, names=names, classes=classes,
                           columns=columns, ...)
          })

## HACK:  seem to need this method
setReplaceMethod( "[",
  signature( x = "timeEvent", i="ANY", j="ANY", value = "timeEvent" ),
  function( x, i, j, ..., value )
  {
    # subscript replace for timeEvent when RHS is also a timeEvent
    # calls the subscript replace with RHS a list
    value <- as( value, class(x))
    if( !identical( value@names, x@names ))
      stop( "Cannot replace -- different structures" )
    len = length(value@names)
    for(k in 1:len){
      x@columns[[k]][i] <- value@columns[[k]]
    }
    x
  })
