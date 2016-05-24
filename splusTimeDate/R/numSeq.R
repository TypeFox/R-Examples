"numericSequence" <- 
function(from, to, by, length.)
{
  ## creation function for sequences
  ## no args -> create default object
  if(missing(from) && missing(to) && missing(by) && missing(length.)) 
    return(new("numericSequence"))
  if(missing(from))
    from <- numeric(0)
  if(missing(to))
    to <- numeric(0)
  if(missing(by))
    by <- numeric(0)
  if(missing(length.))
		length. <- numeric(0)
  lens <- c(length(from), length(to), length(by), length(length.))
  if(any(lens > 1) || (sum(lens) != 3))
    stop(paste(
               "Exactly three of from, to, by, length. must be specified",
               "and they must each be of length 1"))
  if(any(is.na(c(from, to, by, length.))))
    stop("NA values not allowed in sequences")
  ret <- new("numericSequence", from = as(from, "numeric"),
             to = as( to, "numeric" ), by = as(by, "numeric"),
             length = as(length., "integer"))
  ret
}

setAs( "numericSequence", "numeric",
      function( from )
      {
	if( !length( from ))
	  return( numeric(0))

	## OK, must have exactly one of the slot vectors missing, or
	## else we will ignore the length
	## call seq appropriately

	if( length( from@from ) == 0)
	  seq( to = from@to, by = from@by, length = from@length )
	else if( length( from@to ) == 0)
	  seq( from = from@from, by = from@by, length = from@length )
	else if( length( from@by ) == 0)
	  seq( from = from@from, to = from@to, length = from@length )
	else 
	{
	  if( length( from@length ) > 0)
	    warning( "Ignoring length slot in sequence" )
	  seq( from = from@from, to = from@to, by = from@by )
	}
      })

setAs( "numericSequence", "integer", 
      function( from ) as( as( from, "numeric" ), "integer" )
      )

setAs( "numericSequence", "character", 
      function( from ) as( as( from, "numeric" ), "character" )
      )

setAs( "numeric", "numericSequence",
      function( from )
      {
	# convert numeric vector to a sequence, if possible 
	len <- length( from )
	if( !len )
	  return( numericSequence() )
	if( len == 1 )
	  return( numericSequence( from, from, length. = 1 ))

	# verify that it's a valid sequence
	diffs <- from[ -1 ] - from[ - length( from ) ]
	tol <- timeDateOptions( "sequence.tol" )[[1]]
	if( any(abs( diffs - diffs[1] ) > tol ))
	  stop( "Numeric vector is not a regularly spaced sequence" )

	## create sequence from
	numericSequence( from = from[1], by = diffs[1], 
			 length. = len )
      })

setAs( "integer", "numericSequence",
      function( from )
      {
        ## coerce to numeric
        as( from+0 , "numericSequence" )
      })

setMethod( "[", signature( x = "numericSequence", i = "ANY" ), 
           function(x, i, ..., drop = TRUE )
           {
	     y <- as( x, "numeric" )
             y[i, drop=drop]
           })

setMethod( "[[", signature( x = "numericSequence", i = "ANY" ), 
           function(x, ... ) {
		y <- as( x, "numeric" )
                y[[i]]
	})

setReplaceMethod( "[", signature( x = "numericSequence", i = "ANY", j = "ANY"  ),
   function( x, i, j, ..., value ) 
   { 
     x <- as( x, "numeric" )
     x[j] <- value
     x
   })

setReplaceMethod( "[[", signature( x = "numericSequence", i = "ANY" ),
   function( x, i, ..., value ) 
   { 
     x <- as( x, "numeric" )
     x[[i]] <- value
     x
   })

setMethod( "length", signature( x = "numericSequence" ),
   function( x )
   {
   # length for sequences
     slots.missing <- c(length( x@from )==0, 
			length( x@to )==0,
			length( x@by )==0, 
			length( x@length )==0)
     how.many.there <- sum( c(1,1,1,1)[!slots.missing] )

     # if it has a length slot, and missing one of others, that's the length
     if( !slots.missing[4] && ( how.many.there < 4 ))
       return( x@length )
     # if the by slot is 0, then from must be
     # equal to to and it's a length 1 sequence
     # assuming it's valid!
     if( !x@by )
       return( 1 )
     # calculate the length from to/from/by
     as.integer( floor( 1 + (( x@to - x@from ) / x@by )))
   })

setMethod("format", signature="numericSequence",
          function (x, ...) as(x, "character"))


setMethod( "show", "numericSequence", function( object )
{
  if( length( object@from ) )
    cat( "from:  ", object@from, "\n" )
  if( length( object@to ) )
    cat( "to:    ", object@to, "\n" )
  if( length( object@by ) )
    cat( "by:    ", object@by, "\n" )
  if( length( object@length ) )
    cat( "length:", object@length, "\n" )
  nums <- format(as(object,"numeric"))
  if( length(nums) > 5 )
    nums <- c( nums[1:3], "...", nums[length(nums)])
  invisible(print(nums, quote=FALSE))
})

setMethod( "summary", "numericSequence", function( object, ... )
{
  sumry <- numeric(0)
  nms <- character(0)

  from <- object@from
  to <- object@to
  by <- object@by
  len <- object@length

  if( length( from ))
    from <- to - ( len - 1 ) * by
  if( length( to ))
    to <- from + ( len - 1 ) * by 
  if( length( by ))    
    by <- ( to - from ) / ( len - 1 )
  if( length( len ))
    len <- floor( ( to - from ) / by ) + 1 

  sumry <- matrix( c( from, to, by, len), nrow = 1, 
		   dimnames = list( "", c( "From", "To", "By", "Length" )))
  oldClass( sumry ) <- "table"
  sumry
})

setMethod( "Math", "numericSequence", 
	  function( x ) callGeneric( as( x, "numeric" )))

setMethod( "Math2", signature( x = "numericSequence"),  
	  function( x, digits ) callGeneric( as( x, "numeric" ), digits ), 
          )

setMethod( "Summary", signature( x = "numericSequence" ),
	  function(x, ..., na.rm = FALSE) 
            callGeneric( as( x, "numeric" ), ..., na.rm = na.rm ), 
          )
setMethod( "logb", signature( x = "numericSequence" ),
          function(x, base = exp(1))
            callGeneric( as(x, "numeric" ), base = base )
          )
setMethod( "rev", signature( x = "numericSequence" ),
          function(x)
            callGeneric( as(x, "numeric" ) )
          )
setMethod( "Ops", signature( e1 = "numericSequence" ),
	  function( e1, e2 = NULL )
	    callGeneric( as( e1, "numeric" ), e2 ))
setMethod( "Ops", signature( e2 = "numericSequence" ),
	  function( e1, e2 = NULL )
	    callGeneric( e1, as( e2, "numeric" )))
setMethod( "Ops", signature( e1 = "numericSequence", e2 = "numericSequence" ),
	  function( e1, e2 = NULL )
	    callGeneric( as(e1, "numeric"), as( e2, "numeric" )))

setMethod( "is.na", "numericSequence", 
            function( x ) 
            {
	      # valid sequence objects cannot have NA
	      rep(FALSE,length(x))
	    })

setMethod( "match", signature( x = "numericSequence" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	     match( as( x, "numeric" ), table, nomatch, incomparables ), 
           )
setMethod( "match", signature( table = "numericSequence" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	     match( x, as( table, "numeric" ), nomatch, incomparables ), 
           )

setMethod( "unique", signature( x = "numericSequence" ),
           function( x, ... )  x)

setMethod( "duplicated", signature( x = "numericSequence" ),
           function( x, incomparables = FALSE  ) rep( FALSE, length( x )), 
           )

setMethod( "mean", signature( x = "numericSequence" ),
           function(x, trim = 0., na.rm = FALSE, weights = NULL){
	     x <- as( x, "numeric" )
	     callGeneric()
	   })
setMethod( "median", signature( x = "numericSequence" ),
           function(x, na.rm = FALSE){
	     x <- as( x, "numeric" )
	     callGeneric()
	   })
setMethod( "quantile", signature( x = "numericSequence" ),
           function(x, probs = 0:4/4, na.rm = FALSE, ...){
	    x <- as( x, "numeric" )
	    callGeneric()
	  })
setMethod( "sort", signature( x = "numericSequence" ),
           function( x, decreasing = FALSE, ...)
           {
	     # see if it's in reverse order or forward order
	     if( length( x@by ))
	     {
	       if( x@to > x@from )
		 y = x
	       else
	       {
		 tmp <- x@from
		 x@from <- x@to
		 x@to <- tmp
		 y = x
	       }
	     } else {
	       if( x@by >= 0 )
		 y = x
	       else
	       {
		 tmp <- x@from
		 x@from <- x@to
		 x@to <- tmp
		 x@by <- -x@by
		 y = x
	       }
	     }
             if(decreasing) y = rev(y)
             return(y)
	   })

setMethod( "sort.list", signature( x = "numericSequence" ),
           function(x, partial = NULL,  na.last = TRUE, decreasing = FALSE,
                    method = c("shell", "quick", "radix"))
           {
	     ## see if it's in reverse order or forward order
	     if( length( x@by ))
	     {
	       if( x@to > x@from )
		 y = seq(along=x)
	       else
		 y = rev( seq( along=x))
	     } else {
	       if( x@by >= 0 )
		 y = seq( along = x )
	       else
		 y = rev( seq( along=x))
	     }
             if(decreasing) y = rev(y)
             y
	   })


setMethod( "length<-", "numericSequence",
         function( x, value )
       {
	x <- as( x, "numeric" )
	callGeneric()
       })

setMethod( "c", signature(x="numericSequence"),
  function(x, ...){
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
      x <- c(as(x, "numeric"), as( arglist[[1]], "numeric"))
    }
    x
  } 
  )

setMethod("rep", signature(x="numericSequence"),
	  function(x, ...)
  {
    rep(as(x, "numeric"), ...)
  })
   
setMethod("shiftPositions", signature(x="numericSequence"),
	  function(x, k=1)
	  {
	    if(k != round(k)) {
              k <- round(k)
              warning("k is not an integer")
	    }
	    xfrom <- x@from
	    xto <- x@to
        xby <- x@by
        xlen <- x@length
        if(!length(xfrom))
	        xfrom <- xto - (xlen - 1) * xby
	    if(!length(xto))
            xto <- xfrom + (xlen - 1) * xby
        if(!length(xby))
	        xby <- (xto - xfrom) / (xlen - 1)
        if(!length(xlen))
            xlen <- floor((xto - xfrom) / xby) + 1
	    numericSequence(from=xfrom + k*xby, by=xby, length.=xlen)
	  })

setMethod( "diff", signature( x = "numericSequence" ),
    function(x, lag = 1L, differences = 1L, ...) {
        x <- as( x, "numeric" )
        callGeneric()
    })
