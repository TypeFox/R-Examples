setMethod( "-", signature( e1 = "positionsCalendar",
			   e2 = "positionsCalendar" ),
	  function( e1, e2 = NULL )
	  {
	     out <- .Call( "time_time_add", as(e1, "timeDate"), as( e2, "timeDate"),
             -1., "timeSpan")
	     out@format <- as( timeDateOptions("tspan.out.format")[[1]], "character" )
	     out
	   })

setMethod( "-", signature( e1 = "positionsCalendar", e2 = "timeSpan" ),
	  function( e1, e2 = NULL )
	  {
           ret <- .Call( "time_time_add", as(e1, "timeDate"), e2,
               -1., "timeDate")
	       ret@time.zone <- e1@time.zone
	       ret@format <- e1@format
	       ret
	   })

setMethod( "-", signature( e1 = "timeSpan", e2 = "timeSpan" ),
	  function( e1, e2 = NULL )
	  {
            ret <- .Call( "time_time_add", e1, e2, -1., "timeSpan")
	     ret@format <- e1@format
	     ret
	   })

setMethod( "+", signature( e1 = "positionsCalendar", e2 = "timeSpan" ),
	  function( e1, e2 = NULL )
	  {
          ret <- .Call( "time_time_add", as(e1, "timeDate"), e2,
              1., "timeDate")
	      ret@time.zone <- e1@time.zone
	      ret@format <- e1@format
	      ret
	   })

setMethod( "+", signature( e1 = "timeSpan", e2 = "positionsCalendar" ),
	  function( e1, e2 = NULL )
	  {
	     ret <- .Call( "time_time_add", e1, as(e2, "timeDate"), 1., "timeDate")
	     ret@time.zone <- e2@time.zone
	     ret@format <- e2@format
	     ret
	   })


setMethod( "+", signature( e1 = "timeSpan", e2 = "timeSpan" ),
	  function( e1, e2 = NULL )
	  {
	     ret <- .Call( "time_time_add", e1, e2, 1., "timeSpan")
	     ret@format <- e1@format
	     ret
	   })

setMethod( "/", signature( e1 = "timeSpan", e2 = "timeSpan" ),
	  function( e1, e2 = NULL )
	    as( e1, "numeric" ) / as( e2, "numeric" ))


setMethod( "Arith", signature( e1 = "timeSpan", e2 = "numeric" ),
	  function( e1, e2 = NULL )
	  {
	    if(( .Generic != "+" ) && ( .Generic != "-" ))
	      return( callGeneric( as(e1, "numeric"), e2 ))
	    ret <- .Call( "time_num_op", e1, e2, .Generic)
	    ret@format <- e1@format
	    ret
	  })

setMethod( "Arith", signature( e1 = "numeric", e2 = "timeSpan" ),
	  function( e1, e2 = NULL )
	  {
	    if( .Generic == "+" )
	    {
	      ret <- .Call( "time_num_op", e2, e1, .Generic)
	      ret@format <- e2@format
	      ret
	    }
	    else
	      callGeneric( e1, as(e2, "numeric") )
	  })


setMethod( "+", signature( e1 = "timeSpan", e2 = "missing" ),
	   function( e1, e2 = NULL ) e1)

setMethod( "-", signature( e1 = "timeSpan", e2 = "missing" ),
	   function( e1, e2 = NULL )
	   {
	    ret <- .Call( "time_num_op", e1, -1., "*")
	    ret@format <- e1@format
	    ret
	  })

setMethod( "+", signature( e1 = "positionsCalendar", e2 = "numeric" ),
	  function( e1, e2 = NULL )
	  {
	    ret <- .Call( "time_num_op", as(e1, "timeDate"), e2, "+")
	    ret@format <- e1@format
	    ret@time.zone <- e1@time.zone
	    ret
	  })

setMethod( "-", signature( e1 = "positionsCalendar", e2 = "numeric" ),
	  function( e1, e2 = NULL )
	  {
	    ret <- .Call( "time_num_op", as(e1, "timeDate"), e2, "-")
	    ret@format <- e1@format
	    ret@time.zone <- e1@time.zone
	    ret
	  })

setMethod( "+", signature( e1 = "numeric", e2 = "positionsCalendar" ),
	  function( e1, e2 = NULL )
	  {
	    ret <- .Call( "time_num_op", as(e2, "timeDate"), e1, "+")
	    ret@format <- e2@format
	    ret@time.zone <- e2@time.zone
	    ret
	  })

setMethod( "Compare", signature( e1 = "timeSpan", e2 = "numeric" ),
	  function( e1, e2 )
	    callGeneric( as( e1, "numeric" ), e2)
           )

setMethod( "Compare", signature( e1 = "numeric", e2 = "timeSpan" ),
	  function( e1, e2 )
	    callGeneric( e1, as( e2, "numeric" ))
           )

setMethod( "min", signature( x = "positionsCalendar" ),
	   function( x, ..., na.rm = FALSE)
	   {
	     # ignore other arguments -- generic takes care of them
	     ret <- .Call( "time_range", as(x, "timeDate"),  na.rm)
	     ret@format <- x@format
	     ret@time.zone <- x@time.zone
	     ret[1]
	   })

setMethod( "max", signature( x = "positionsCalendar" ),
	   function( x, ..., na.rm = FALSE)
	   {
	     # ignore other arguments -- generic takes care of them
	     ret <- .Call( "time_range", as(x, "timeDate"), na.rm)
	     ret@format <- x@format
	     ret@time.zone <- x@time.zone
	     ret[2]
	   })

setMethod( "range", signature( x = "positionsCalendar" ),
	   function( x, ..., na.rm = FALSE)
	   {
	     # ignore other arguments -- generic takes care of them
	     ret <- .Call( "time_range", as(x, "timeDate"), na.rm)
	     ret@format <- x@format
	     ret@time.zone <- x@time.zone
	     ret
	   })

setMethod( "min", signature( x = "timeSpan" ),
	   function( x, ..., na.rm = FALSE)
	   {
         z <- list(...)
         for(i in z) x <- c(x, i)
	     # ignore other arguments -- generic takes care of them
	     ret <- .Call( "time_range", x, na.rm)
	     ret@format <- x@format
	     ret[1]
	   })

setMethod( "max", signature( x = "timeSpan" ),
	   function( x, ..., na.rm = FALSE)
	   {
         z <- list(...)
         for(i in z) x <- c(x, i)
	     # ignore other arguments -- generic takes care of them
	     ret <- .Call( "time_range", x, na.rm)
	     ret@format <- x@format
	     ret[2]
	   })

setMethod( "range", signature( x = "timeSpan" ),
	   function( x, ..., na.rm = FALSE)
	   {
         z <- list(...)
         for(i in z) x <- c(x, i)
	     # ignore other arguments -- generic takes care of them
	     ret <- .Call( "time_range", x, na.rm)
	     ret@format <- x@format
	     ret
	   })


setMethod( "sum", signature( x = "timeSpan" ),
	   function( x, ..., na.rm = FALSE)
	   {
         z <- list(...)
         for(i in z) x <- c(x, i)
	     # ignore other arguments -- generic takes care of them
	     ret <- .Call( "time_sum", x, na.rm, FALSE)
	     ret@format <- x@format
	     ret
	   })

setMethod( "cumsum", signature( x = "timeSpan" ),
	   function( x )
	   {
	     ret <- .Call( "time_sum", x, FALSE, TRUE)
             ret@format <- x@format
	     ret
	   })

setMethod( "mean", signature( x = "positionsCalendar" ),
	   function(x, trim = 0.0, na.rm = FALSE, weights = NULL)
	   {
	     ret <- mean( as( x, "numeric" ), trim, na.rm, weights)
	     ret <- as( ret, "timeDate" )
	     ret@format <- x@format
	     ret@time.zone <- x@time.zone
	     ret
	   })

setMethod( "mean", signature( x = "timeSpan" ),
	   function(x, trim = 0.0, na.rm = FALSE, weights = NULL)
	   {
	     ret <- mean( as( x, "numeric" ), trim, na.rm, weights )
	     ret <- as( ret, "timeSpan" )
	     ret@format <- x@format
	     ret
	   })

setMethod( "median", signature( x = "positionsCalendar" ),
	   function(x, na.rm = FALSE)
	   {
	     ret <- median( as( x, "numeric" ), na.rm)
	     ret <- as( ret, "timeDate" )
	     ret@format <- x@format
	     ret@time.zone <- x@time.zone
	     ret
	   })

setMethod( "median", signature( x = "timeSpan" ),
	   function(x, na.rm = FALSE)
	   {
	     ret <- median( as( x, "numeric" ), na.rm )
	     ret <- as( ret, "timeSpan" )
	     ret@format <- x@format
	     ret
	   })


setMethod( "quantile", signature( x = "positionsCalendar" ),
	   function(x, probs = 0:4/4, na.rm = FALSE, ...)
	   {
	     ret <- quantile( as( x, "numeric" ), probs, na.rm, ...)
	     ret <- as( ret, "timeDate" )
	     ret@format <- x@format
	     ret@time.zone <- x@time.zone
	     ret
	   })

setMethod( "quantile", signature( x = "timeSpan" ),
	   function(x, probs = 0:4/4, na.rm = FALSE, ...)
	   {
	     ret <- quantile( as( x, "numeric" ), probs, na.rm, ...)
	     ret <- as( ret, "timeSpan" )
	     ret@format <- x@format
	     ret
	   })


setMethod("cor", signature( x = "positionsCalendar", y = "ANY" ),
	  function(x, y = NULL, use = "everything",
                   method = c("pearson", "kendall", "spearman"))
	  {
	    x <- as( x, "numeric" )
	    callGeneric()
	  })

setMethod("cor", signature( x = "ANY", y = "positionsCalendar" ),
	  function(x, y = NULL, use = "everything",
                   method = c("pearson", "kendall", "spearman"))
	  {
	    y <- as( y, "numeric" )
	    callGeneric()
	  })

setMethod("cor", signature( x = "positionsCalendar", y = "positionsCalendar" ),
	  function(x, y = NULL, use = "everything",
                   method = c("pearson", "kendall", "spearman"))
	  {
	    x <- as( x, "numeric" )
	    y <- as( y, "numeric" )
	    callGeneric()
	  })

setMethod("cor", signature( x = "timeSpan", y="timeSpan" ),
	  function(x, y = NULL, use = "everything",
                   method = c("pearson", "kendall", "spearman"))
	  {
	    x <- as( x, "numeric" )
	    y <- as( y, "numeric" )
	    callGeneric()
	  })

setMethod("cor", signature( x = "timeSpan", y="ANY" ),
	  function(x, y = NULL, use = "everything",
                   method = c("pearson", "kendall", "spearman"))
	  {
	    x <- as( x, "numeric" )
	    callGeneric()
	  })

setMethod("cor", signature( x = "ANY", y="timeSpan" ),
	  function(x, y = NULL, use = "everything",
                   method = c("pearson", "kendall", "spearman"))
	  {
	    y <- as( y, "numeric" )
	    callGeneric()
	  })


setMethod("var", signature( x = "positionsCalendar" ),
	  function(x, y, na.rm = FALSE, use) 
	  {
	    x <- as( x, "numeric" )
	    callGeneric()
	  })


setMethod("var", signature( x = "timeSpan" ),
	  function(x, y, na.rm = FALSE, use) 
	  {
	    x <- as( x, "numeric" )
	    callGeneric()
	  })


setMethod("var", signature( y = "positionsCalendar" ),
	  function(x, y, na.rm = FALSE, use) 
	  {
	    y <- as( y, "numeric" )
	    callGeneric()
	  })


setMethod("var", signature( y = "timeSpan" ),
	  function(x, y, na.rm = FALSE, use) 
	  {
	    y <- as( y, "numeric" )
	    callGeneric()
	  })

setMethod( "diff", signature( x = "positionsCalendar" ),
	   function( x, lag = 1, differences = 1 )
	   {
	     xlen <- length( x )
	     lag <- round( lag )
	     differences <- round( differences )
	     if(lag < 1 | differences < 1 )
	       stop("Bad value for lag or differences")
	     if(lag * differences >= xlen)
	       return(x[0])
	     s <- 1:lag
	     for(i in 1:differences)
	       x <- x[ - s] - x[ - (length(x) + 1 - s)]
	     x
	   })

setMethod( "diff", signature( x = "timeSpan" ),
	   function( x, lag = 1, differences = 1 )
	   {
	     xlen <- length( x )
	     lag <- round( lag )
	     differences <- round( differences )
	     if(lag < 1 | differences < 1 )
	       stop("Bad value for lag or differences")
	     if(lag * differences >= xlen)
	       return(x[0])
	     s <- 1:lag
	     for(i in 1:differences)
	       x <- x[ - s] - x[ - (length(x) + 1 - s)]
	     x
	   })



setMethod( "match", signature( x = "positionsCalendar", table = "positionsCalendar" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     # use numeric matching
	     if( is( incomparables, "positionsCalendar" ))
	       incomparables <- as( incomparables, "numeric" )
	     match( as( x, "numeric" ), as( table, "numeric" ),
		    nomatch, incomparables )
	   })

setMethod( "match", signature( x = "character", table = "positionsCalendar" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	     match( as( x, "timeDate" ), table, nomatch, incomparables ),
	  )

setMethod( "match", signature( x = "positionsCalendar", table = "character" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     if( is( incomparables, "positionsCalendar" ))
	       incomparables <- as( incomparables, "character" )
	     match( as( x, "character" ), table, nomatch, incomparables )
	   })

setMethod( "match", signature( x = "positionsCalendar" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     warning( paste( "Cannot match time with", class( table ),
			  "result will be nomatch" ))
	     rep( nomatch, length( x ))
	   })

setMethod( "match", signature( table = "positionsCalendar" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     warning( paste( "Cannot match", class( x ), "with time",
			  "result will be nomatch" ))
	     rep( nomatch, length( x ))
	   })

setMethod( "match", signature( x = "timeSpan", table = "timeSpan" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     # use numeric matching
	     if( is( incomparables, "timeSpan" ))
	       incomparables <- as( incomparables, "numeric" )
	     match( as( x, "numeric" ), as( table, "numeric" ),
		    nomatch, incomparables )
	   })

setMethod( "match", signature( x = "character", table = "timeSpan" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	     match( as( x, "timeSpan" ), table, nomatch, incomparables )
	   )

setMethod( "match", signature( x = "timeSpan", table = "character" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     if( is( incomparables, "timeSpan" ))
	       incomparables <- as( incomparables, "character" )
	     match( as( x, "character" ), table, nomatch, incomparables )
	   })

setMethod( "match", signature( x = "timeSpan" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     warning( paste( "Cannot match timeSpan with", class( table ),
			  "result will be nomatch" ))
	     rep( nomatch, length( x ))
	   })

setMethod( "match", signature( table = "timeSpan" ),
           function(x, table, nomatch = NA, incomparables = FALSE)
	   {
	     warning( paste( "Cannot match", class( x ), "with timeSpan",
			  "result will be nomatch" ))
	     rep( nomatch, length( x ))
	   })


##setMethod( "ordered", signature( x = "positionsCalendar" ),
##	  function(x, levels = sort(unique(x)),
##		   labels = as(levels,"character"), exclude = NA)
##	  {
##	    y <- factor( x, levels, labels, exclude )
##	    oldClass( y ) <- c( "ordered", "factor" )
##	    y
##	  })


##setMethod( "ordered", signature( x = "timeSpan" ),
##	  function(x, levels = sort(unique(x)),
##		   labels = as(levels,"character"), exclude = NA)
##	  {
##	    y <- factor( x, levels, labels, exclude )
##	    oldClass( y ) <- c( "ordered", "factor" )
##	    y
##	  })

setMethod( "cut", signature( x = "positionsCalendar" ),
	   function(x, breaks, labels, include.lowest = FALSE, factor.result = FALSE, right=FALSE)
	  {
	    if( is( breaks, "character" ) && ( length(breaks) == 1 ))
	    {
	      # breaks must be something you could use as by
	      # argument for sequence
	      if( length(x)) {
            bybrk <- breaks
            rng <- range(x, na.rm=TRUE)
            breaks <- timeSeq(from=rng[1], to=rng[2], by=bybrk,
                              align.by=TRUE, extend=TRUE, format=x@format,
                              zone=x@time.zone)
            if( !include.lowest ) {
                if( !right && ( breaks[length(breaks)] == rng[2] ))
                  breaks <- c( breaks,
                                   breaks[length(breaks)] + timeRelative(by=bybrk))
                  else if( right && ( breaks[1] == rng[1] ))
                    breaks <- c(breaks[1] - timeRelative(by=bybrk),
                                     breaks)
              }
          } else {
            breaks <- timeDate()
          }
	    } else if( is(breaks, "numeric") && length( breaks) == 1 ) {
        if(breaks < 1)
          stop("Must specify at least one interval")
        if(!missing(labels) && length(labels) != breaks)
          stop("Number of labels must equal number of intervals")

	      # cut into N equal intervals
        if( length(x)) {
          rng <- range(x, na.rm=TRUE)
          rng[is.na(rng)] <- timeDate(julian=0)
          if( rng[1] == rng[2] ) 
            rng[2] <- rng[2] + timeSpan("1MS" )
          breaks <- timeSeq(from=rng[1], to=rng[2],
                            length.out=breaks + 1,
                            format=x@format, zone=x@time.zone)
          if( !include.lowest ) {
            if( !right )
              breaks[length(breaks)] <-
                breaks[length(breaks)] + timeSpan("1MS")
              else
                breaks[1] <- breaks[1] - timeSpan("1MS")
          }
        } else {
          breaks <- timeDate()
        }
      } else {
        breaks <- as( breaks, "timeDate" )
      }
	    if(!right)
	      labpaste <- c(" thru ", "-")
	      else labpaste <- c("+ thru ", "")
	    if(missing(labels))
	      labels <- paste(as(breaks[ - length(breaks)], "character"),
                        labpaste[1], as(breaks[-1], "character"),
                        labpaste[2], sep = "")
	    cut( as(x,"numeric"), breaks = as( breaks, "numeric" ),
          labels = labels, include.lowest = include.lowest, right = right,
          factor.result=factor.result )
	  })

setMethod( "cut", signature( x = "timeSpan" ),
	   function(x, breaks, labels, include.lowest = FALSE, factor.result = FALSE, right = FALSE)
	  {
      if( is(breaks, "numeric") && length( breaks) == 1 ) {
        if(breaks < 1)
          stop("Must specify at least one interval")
        if(!missing(labels) && length(labels) != breaks)
          stop("Number of labels must equal number of intervals")

        # cut into N equal intervals
        include.lowest <- TRUE
        right = TRUE
        if(length(x)) {
          rng <- range(x, na.rm=TRUE)
          rng[is.na(rng)] <- timeSpan(julian=0)
          if( rng[1] == rng[2] ) 
            rng[2] <- rng[2] + timeSpan("1MS" )
          rng <- as.double(rng)
          breaks <- timeSpan(julian=seq(from=rng[1], to=rng[2],
            length.out=breaks + 1))
          if( !include.lowest ) {
            if( !right )
              breaks[length(breaks)] <-
                breaks[length(breaks)] + timeSpan("1MS")
              else
                breaks[1] <- breaks[1] - timeSpan("1MS")
          }
        } else {
          breaks <- timeSpan()
        }
      } else { 
        breaks <- as( breaks, "timeSpan" )
      }
	    if(!right)
	      labpaste <- c(" thru ", "-")
	      else labpaste <- c("+ thru ", "")
	    if(missing(labels))
	      labels <- paste(as(breaks[ - length(breaks)], "character"),
                        labpaste[1], as(breaks[-1],  "character"),
                        labpaste[2], sep = "")
	    cut( as(x,"numeric"), breaks = as( breaks, "numeric" ),
          labels = labels, include.lowest = include.lowest, right = right,
          factor.result = factor.result )
	  })

