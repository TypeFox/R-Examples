"timeRelative" <- 
function(x, holidays., by, k.by = 1, align.by = FALSE, week.day = NULL)
{
  obj <- new("timeRelative")
  if(!missing(x)) {
    x <- as(x, "character")
    obj@Data <- x
    if(!missing(by) || !missing(k.by) || !missing(align.by))
      warning("by, k.by, and align.by arguments ignored")
  }
  else if(!missing(by)) {
    k.by <- as(k.by, "integer")
    if(k.by == 0)
      stop("k.by cannot be zero")
    direction = 1
    if(k.by < 0) {
      direction = -1
      k.by <- abs(k.by)
    }
    if( by == "wk" ) by <- "weeks"
    switch(by,
           milliseconds = str <- "ms",
           ms = str <- "ms",
           seconds = str <- "sec",
           sec = str <- "sec",
           minutes = str <- "min",
           min = str <- "min",
           hours = str <- "hr",
           hr = str <- "hr",
           days = str <- "day",
           day = str <- "day",
           weekdays = str <- "wkd",
           wkd = str <- "wkd",
           bizdays = str <- "biz",
           biz = str <- "biz",
           weeks = {
             if(is.null(week.day)) {
               if(align.by) {
                 str <- "day"
                 k.by <- k.by * 7
               }
               else str <- "wk"
             }
             else {
               if(is(week.day, "numeric"))
                 move.to.day <- week.day + 1
               else move.to.day <-
                 charmatch(tolower( week.day),
                           tolower( timeDateOptions("time.day.name"
                                             )[[1]]), nomatch = 0)
               if(length(move.to.day) > 1) {
                 move.to.day <- move.to.day[
                                            1]
                 warning("Only first week.day string used"
                         )
               }
               if((move.to.day < 1) || (move.to.day >
                                        7))
                 stop(paste(
                            "Do not know how to align to day",
                            week.day))
               str <- c("sun", "mon", "tue", "wed",
                        "thu", "fri", "sat")[
                                             move.to.day]
             }
           }
           ,
           months = str <- "mth",
           mth = str <- "mth",
           quarters = str <- "qtr",
           qtr = str <- "qtr",
           years = str <- "yr",
           yr = str <- "yr",
           sun = str <- "sun",
           mon = str <- "mon",
           tue = str <- "tue",
           wed = str <- "wed",
           thu = str <- "thu",
           fri = str <- "fri",
           sat = str <- "sat",
           stop("Unknown by argument"))
    str <- paste(if(direction > 0) "+" else "-", if(align.by) "a"
    else "", k.by, str, sep = "")
    obj@Data <- str
  }
  if(!missing(holidays.)) {
    if(!is(holidays., "positionsCalendar"))
      holidays. <- as(holidays., "timeDate")
    obj@holidays <- holidays.
  } else obj@holidays = timeDate()
  obj
}

setMethod( "initialize", "timeRelative",
          function( .Object, Data, holidays, ...){
            if(missing(Data))
              Data <- character(0)
            if(missing(holidays))
              holidays <- new("timeDate")
            callNextMethod(.Object, Data=Data, holidays=holidays, ...)
          })

setAs( "character", "timeRelative", 
       function( from ) timeRelative( from )
       )
setAs( "timeRelative", "character", 
       function( from ) from@Data
       )

setMethod( "format", signature(x="timeDate"),
          function(x, ...)
          {
            as(x, "character")
          })

setMethod( "show", "timeRelative", 
	   function( object )
	   {
	     if(length(object))
	       show( object@Data )
	     else
	       cat("timeRelative()\n")
	   } )

setMethod( "summary", "timeRelative", function( object, ... )
	  {
	    ret <- matrix( c( length( object@Data ), "timeRelative" ),
			   nrow=1,
			   dimnames = list( "", 
			     c( "Length", "Class" )))
	    oldClass( ret ) <- c( "table", "matrix" )
	    ret
	  }  )

setMethod( "[", signature(x = "timeRelative", i = "ANY"),
	  function( x, i, ..., drop = TRUE ) 
	  { 
	    x@Data <- x@Data[i, drop=drop]
	    x
	  } )

setMethod( "[[", signature(x = "timeRelative", i = "ANY"),
	  function( x, i, ... ) 
	  { 
	    x@Data <- x@Data[[i]]
	    x
	  } )

setReplaceMethod( "[", signature( x = "timeRelative", i="ANY", j="ANY"),
	  function( x, i, j, ..., value ) 
	  { 
	    tmp <- x@Data
	    tmp[i] <- as( value, "character" )
	    x@Data <- tmp
	    x
	  } )

setReplaceMethod( "[[", signature( x = "timeRelative" , i = "ANY"),
	  function( x, ..., value ) 
	  { 
	    tmp <- x@Data
	    tmp[[...]] <- as( value, "character" )
	    x@Data <- tmp
	    x
	  } )

setMethod( "length", "timeRelative", function(x) length( x@Data ),
	  )

setReplaceMethod( "length", "timeRelative", 
		 function(x,value)
		 {
		   length(x@Data) <- value
		   x
		 } )

setMethod( "is.na", "timeRelative", function(x) is.na( x@Data ),
	  )


setMethod( "-", signature( e1 = "timeRelative", e2 = "missing" ),
	   function( e1, e2 )
	   {
	     e1@Data <- timerel.chsigns( e1@Data )
	     e1
	   } )

setMethod( "+", signature( e1 = "timeRelative", e2 = "timeRelative" ),
	   function( e1, e2 ) 
	   { 
	     e1@Data <- paste( e1@Data, e2@Data )
	     e1 
	   } )

setMethod( "-", signature( e1 = "timeRelative", e2 = "timeRelative" ),
	   function( e1, e2 ) 
	   { 
	     e1@Data <- paste( e1@Data, timerel.chsigns(e2@Data))
	     e1 
	   } )

setMethod( "*", signature( e1 = "timeRelative", e2 = "numeric" ),
	   function( e1, e2 )
	   {
             if (any(e2 != as.integer(e2), na.rm = TRUE))
                 stop("Can only multiply a timeRelative object by an integer")
	     e1@Data <- sapply( e1@Data, 
			   function(x, n ) paste( rep( x, n ), collapse=" "),
			   as.integer(e2) )
	     e1 
	   } )

setMethod( "*", signature( e1 = "numeric", e2 = "timeRelative" ),
	   function( e1, e2 ) e2 * e1  )

setMethod( "Ops", signature( e1 = "timeRelative" ),
	   function( e1, e2 ) 
	     stop( "Function does not make sense for class timeRelative" )
	   )
setMethod( "Ops", signature( e2 = "timeRelative" ),
	   function( e1, e2 ) 
	     stop( "Function does not make sense for class timeRelative" )
	   )

setMethod( "Math", "timeRelative", 
	  function( x ) stop( "Function does not make sense for class timeRelative" )
	  )

setMethod( "logb", signature( x = "timeRelative" ),
	  function( x, base ) stop( "Function does not make sense for class timeRelative" )
	  )

setMethod( "Math2", signature( x = "timeRelative"),  
	  function( x, digits ) stop( "Function does not make sense for class timeRelative" )
	  )

setMethod( "Summary", signature( x = "timeRelative" ),
	  function( x, ..., na.rm = FALSE ) stop( "Function does not make sense for class timeRelative" )
	  )

setMethod( "+", signature( e1 = "positionsCalendar", e2 = "timeRelative" ),
	   function( e1, e2 ) 
	   {
	     e1 <- as( e1, "timeDate" )
	     tmp <- .Call( "time_rel_add",
                          as(e1, "timeDate"),
                          as(e2@Data, "character"),
                          as(sort(e2@holidays), "timeDate"),
                          timeZoneList())
	     e1@columns <- tmp@columns
	     e1
	   } ) 

setMethod( "+", signature( e1 = "timeRelative", e2 = "positionsCalendar" ),
	   function( e1, e2 ) e2 + e1  )

setMethod( "-", signature( e1 = "positionsCalendar", e2 = "timeRelative" ),
	   function( e1, e2 ){ 
	     e2@Data <- timerel.chsigns( e2@Data ) 
	     e1 + e2
	   } )


setMethod( "match", signature( x = "timeRelative" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	      match( x@Data, table, nomatch, incomparables )
	   )
setMethod( "match", signature( table = "timeRelative" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	      match( x, table@Data, nomatch, incomparables )
	   )

setMethod( "c", signature(x="timeRelative"),
  function(x, ...){
    ## Concatenates the data slot of the two timeRelative objects, if
    ## they are of the same class or can be coerced to be the same.
    arglist <- list(...)
    if(length(arglist)==0) return(x)
    lens <- sapply(arglist, length)
    if(!any(lens > 0)) {
      return(x)
    }
    arglist <- arglist[lens > 0]
    if(length(arglist) > 1)
      c(c(x, arglist[[1]]), do.call("c", arglist[-1]))
    else {
      y <- as( arglist[[1]], class( x ))
      x@Data <- c(x@Data, y@Data)
      x
    } 
  })

setMethod( "duplicated", "timeRelative", 
	  function( x, incomparables=FALSE)
	   duplicated( x@Data, incomparables )
	   )

"timerel.chsigns" <- 
function(x)
{
  ## Undocumented function to change the sign for relative time class.
  ## Basically just changes + to - and - to + in a character vector.
  maxn <- max(nchar(x))
  charmat <- sapply(x, "substring", 1:maxn, 1:maxn)
  plus <- charmat == "+"
  minus <- charmat == "-"
  charmat[plus] <- "-"
  charmat[minus] <- "+"
  apply(charmat, 2, "paste", collapse = "")
}

"is.monthend" <- 
function(x)
{
	nextday <- x + timeRelative(by = "days", k.by = 1)
	mdy(nextday)$day == 1
}
