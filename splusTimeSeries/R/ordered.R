
setAs( "seriesVirtual", "numeric",
      function(from) as(from@data, "numeric"))
setAs( "seriesVirtual", "complex",
      function(from) as(from@data, "complex"))
setAs( "seriesVirtual", "character",
      function(from) as(from@data, "character"))
setAs( "seriesVirtual", "matrix",
      function(from) as(from@data, "matrix"))
setAs( "seriesVirtual", "logical",
      function(from) as(from@data, "logical"))
setAs( "seriesVirtual", "integer",
      function(from) as(from@data, "integer"))
setAs( "seriesVirtual", "vector",
      function(from) as(from@data, "vector"))

"as.data.frame.timeSeries" <- 
function(x, ...) {
    V1 <- x@data
    ret <- as.data.frame(V1, ...)
    row.names(ret) <- as(x@positions, "character")
    ret
}
"as.data.frame.signalSeries" <-
function(x, ...) {
    ret <- as.data.frame(x@data, ...)
    row.names(ret) <- as(x@positions, "character")
    ret
}
setAs("timeSeries", "data.frame",
      function(from) as.data.frame.timeSeries(from))
setAs("signalSeries", "data.frame",
      function(from) as.data.frame.signalSeries(from))

setMethod( "show", "seriesVirtual", function( object )
{
  # print out the data by making a table
  if( !numRows( object@data ))
  {
    cat( "NULL data\n" )
    return()
  }

  if(is.data.frame( object@data)) {
      cols <- unlist( lapply( 1:numCols(object@data),
          function(i, x ) format(subscript2d( x, , i )), object@data ))
  } else {
      cols <- unlist( lapply( 1:numCols(object@data),
          function(i, x ) subscript2d( x, , i ),
                         object@data ))
  }

  cols <- c( as( object@positions, "character" ), cols )

  cols <- matrix( cols, ncol = numCols( object@data ) + 1 )
  colnames <- names(object@data)
  if( is.null( colnames ) || !length( colnames ))
    colnames <- 1:numCols(object@data)
  dimnames( cols ) <- list( rep("",numRows(cols)),
			   c( "Positions", colnames ))

  oldClass( cols ) <- "table"
  show( cols )
})

setMethod( "summary", "seriesVirtual", function( object, ... )
{
  ps <- summary( object@positions )
  # assume summary() always returns a 1-row table; want matrix
  oldClass(ps) <- NULL
  ds <- lapply( 1:numCols(object@data),
	        function(i, x)
	        {
		  ret <- summary(subscript2d(x,,i))
	          oldClass(ret) <- NULL
	          ret
		}, object@data )

  lens <- sapply( ds, "numCols" )
  lens <- c( numCols(ps), lens )
  maxlen <- max( lens )

  # make all summaries into name : value character strings

  pastenames <- function( tab, newlen )
  {
    ret <- paste( format( dimnames(tab)[[2]] ), ":",  format( tab ),
		  "  ", sep = "" )
    length( ret ) <- newlen
    ret
  }

  allcols <- cbind( pastenames( ps, maxlen ),
		    sapply( ds, pastenames, maxlen ))

  colnames <- names(object@data)
  if( is.null( colnames ))
    colnames <- 1:numCols(object@data)
  dimnames( allcols ) <- list( rep("",maxlen),
			   c( "Positions", colnames ))
  oldClass( allcols ) <- "table"
  allcols
})

setMethod( "nrow", "seriesVirtual", function( x ) numRows( x@data ))
setMethod( "numRows", "seriesVirtual",
function( x ) numRows( x@data ) )

setMethod( "ncol", "seriesVirtual", function( x ) numCols( x@data ) )
setMethod( "numCols", "seriesVirtual", function( x ) numCols( x@data ) )

setMethod("dim", "seriesVirtual", function(x) dim(x@data ) )

setMethod( "length", "seriesVirtual", function( x ) length( x@data ) )

setMethod( "seriesLength", "seriesVirtual", function( x ) length( x@positions ) )

setMethod( "format", "seriesVirtual", function( x, ... ) format( x@data, ... ) )

setMethod("duplicated", "seriesVirtual", function( x, incomparables=FALSE) {
	x <- x@data
	callGeneric()
} )

setMethod("dim", "seriesVirtual", function(x) dim(x@data))

setMethod("row.names", "seriesVirtual", function(x) x@positions )

setMethod("names", "seriesVirtual", function(x) dimnames(x@data)[[2]] )

setMethod("dimnames", "seriesVirtual",
   function(x) list(row.names(x), dimnames(x@data)[[2]]))

setReplaceMethod("row.names", "seriesVirtual",
   function(x, value )
   {
     positions(x) <- value
     x
   })

setReplaceMethod("names", "seriesVirtual",
   function(x, value )
   {
     dimnames(x@data)[[2]] <- value
     x
   })

## taken from R version for data frame
setReplaceMethod("dimnames", "seriesVirtual",
   function(x, value )
   {
     d <- dim(x)
     if(!is.list(value) || length(value) != 2L)
       stop("invalid 'dimnames' given for data frame")
     ## do the coercion first, as might change length
     value[[1L]] <- as.character(value[[1L]])
     value[[2L]] <- as.character(value[[2L]])
     if(d[[1L]] != length(value[[1L]]) || d[[2L]] != length(value[[2L]]))
       stop("invalid 'dimnames' given for data frame")
     row.names(x) <- value[[1L]] # checks validity
     names(x) <- value[[2L]]
     x
   })

setReplaceMethod( "numRows", signature(x="seriesVirtual"),
                  function(x, value )
                  {
		    numRows(x@positions) <- value
		    numRows(x@data ) <- value
		    x
		  })

setReplaceMethod( "numCols", "seriesVirtual",
                  function(x, value )
                  {
		    numCols(x@data ) <- value
		    x
		  })

setMethod("rowIds", "seriesVirtual", function(x) x@positions )

setMethod("colIds", "seriesVirtual", function(x) colIds(x@data) )

setReplaceMethod("rowIds", "seriesVirtual",
   function(x, value )
   {
     positions(x) <- value
     x
   })

setReplaceMethod("colIds", "seriesVirtual",
   function(x, value )
   {
     colIds(x@data) <- value
     x
   })


setMethod( "start", "seriesVirtual",
function( x, ... )
  {
    if( length( x@start.position ))
      x@start.position
    else if( length( x@positions ))
      x@positions[1]
    else
      stop( "No start position" )
  })

setMethod( "end", "seriesVirtual",
function( x, ... )
  {
    if( length( x@end.position ))
      x@end.position
    else if( length( x@positions ))
      x@positions[length( x@positions )]
    else
      stop( "No end position" )
  })

## Hack - this shouldn't be needed by [ doesn't work for sequences
setMethod( "start", "signalSeries",
function( x, ... )
  {
    if( length( x@start.position ))
      x@start.position
    else if( length( x@positions ))
      as(x@positions, "numeric")[1]
    else
      stop( "No start position" )
  })

## Hack - this shouldn't be needed by [ doesn't work for sequences
setMethod( "end", "signalSeries",
function( x, ... )
  {
    if( length( x@end.position ))
      x@end.position
    else if( length( x@positions ))
      as(x@positions, "numeric")[length( x@positions )]
    else
      stop( "No end position" )
  })

setMethod("subscript2d", "seriesVirtual",
	  function(x,i,j)
	  {
	    if( !missing(i) && !missing(j))
	      return(x[i,j,drop=FALSE])
	    if( !missing(i))
	      return(x[i,,drop=FALSE])
	    if( !missing(j))
	      return(x[,j,drop=FALSE])
	    x[,,drop=FALSE]
	  })

setMethod( "[", signature(x="seriesVirtual"),
          function( x, i, j, drop = TRUE )
{
  ## drop argument is ignored when subscripting time series!
  has.i <- !missing(i)
  has.j <- !missing(j)
  if(!has.i && !has.j)
    return(x)

  # if i,j args are time series, convert them
  if( has.i && is( i, "series" ))
    i <- seriesData( i )
  if( has.j && is( j, "series" ))
    j <- seriesData( j )

## handle timeEvent subscripting
  if( has.i && is(i, "timeEvent")) {
	mypos <- x@positions
	nevent <- length(i)
	idx <- logical(length(mypos))
	for (k in 1:nevent){
		idx <- idx | ((i@columns[[1]][k] <= mypos) & (i@columns[[2]][k] >= mypos))
	}
	i <- idx
  }

  ## HACK: this is because [ operator doesn't work
  if( has.i ){
    if(inherits(x@positions, "positionsNumeric"))
      x@positions <- as(x@positions, "numeric")
    x@positions <- x@positions[i]
  }
  oldnames <- names( x )
  ## this next line allows subscripting by names
  names( oldnames ) <- oldnames

  if(has.i && has.j)
    x@data <- asSeriesData( subscript2d(x@data, i, j) )
  else if( !has.i && has.j )
    x@data <- asSeriesData( subscript2d(x@data, , j) )
  else if( has.i && !has.j )
    x@data <- asSeriesData( subscript2d(x@data, i, ) )

  if( has.j )
    newnames <- oldnames[j]
  else
    newnames <- oldnames
  names( newnames ) <- NULL
  nc <- numCols( x@data )
  # may have ended up with a null data matrix here
  if( !is.data.frame(x@data) && !is.null(newnames))
    row.names(x@data) <- character(0)
  if( nc && !is.null(newnames))
    names( x ) <- newnames[1:nc]
  x
})

setReplaceMethod( "[", signature( x = "seriesVirtual", value = "vector" ),
function( x, i, j, ..., value )
{
  has.i <- !missing(i)
  has.j <- !missing(j)

  # if i,j args are time series, convert them
  if( has.i && is( i, "seriesVirtual" ))
    i <- seriesData( i )
  if( has.j && is( j, "seriesVirtual" ))
    j <- seriesData( j )

## handle timeEvent subscripting
  if( has.i && is(i, "timeEvent")) {
	mypos <- x@positions
	nevent <- length(i)
	idx <- matrix(FALSE, nrow=nevent, ncol=length(mypos))
	for (k in 1:nevent){
          idx[k,] <- (i@columns[[1]][k] <= mypos) & (i@columns[[2]][k] >= mypos)
	}
	i <- as.logical(colSums(idx))
  }

  if( has.i && has.j )
    subscript2d(x@data, i, j) <- value
  else if( !has.i && has.j )
    subscript2d(x@data, , j) <- value
  else if( has.i && !has.j )
    subscript2d(x@data, i, ) <- value
  else # both null
    x@data <- value

  x
})

setReplaceMethod( "[", signature( x = "seriesVirtual",
                                 value = "seriesVirtual" ),
function( x, i, j, ..., value )
{
  ## no check is done on the positions of the RHS in this case!
  x[i, j] <- value@data
  x
})

setMethod( "Math", "seriesVirtual",
   function( x )
   {
     x@data <- asSeriesData( callGeneric( x@data ))
     x
   })


setReplaceMethod("subscript2d", "seriesVirtual",
function(x,i,j,value)
  {
    if( !missing(i) && !missing(j))
      x[i,j] <- value
    else if( missing(i) && !missing(j))
      x[,j] <- value
    else if( !missing(i) && missing(j))
      x[i,] <- value
    else
      x[] <- value
    x
  })


setMethod( "Math2", signature( x = "seriesVirtual" ),
   function( x, digits )
   {
     x@data <- asSeriesData( callGeneric( x@data, digits ))
     x
   })

setMethod( "Summary", signature( x = "seriesVirtual" ),
           function( x, ..., na.rm = FALSE )
	     callGeneric( x@data, na.rm = na.rm ))

setMethod("is.na", "seriesVirtual",
  function(x)
  {
    x@data <- is.na(x@data)
    x
  })

setMethod( "Ops", signature( e1 = "seriesVirtual" ),
           function( e1, e2 = NULL )
           {
	     e1@data <- asSeriesData( callGeneric(e1@data, e2))
	     e1
	   })


setMethod( "Ops", signature( e1 = "seriesVirtual", e2="missing" ),
           function( e1, e2 = NULL )
           {
	     e1@data <- asSeriesData( callGeneric(e1@data))
	     e1
	   })

setMethod( "Ops", signature( e2 = "seriesVirtual" ),
           function( e1, e2 = NULL )
           {
	     e2@data <- asSeriesData( callGeneric(e1, e2@data))
	     e2
	   })

setMethod( "Ops", signature( e1 = "seriesVirtual", e2 = "seriesVirtual" ),
           function( e1, e2 = NULL )
           {
	     if( is( positions( e1 ), "positionsCalendar" ) !=
		 is( positions( e2 ), "positionsCalendar" ))
	       stop( "Cannot operate on numeric and calendar series" )

	     # make union of positions
	     ##newpos <- sort(union( positions( e1 ), positions( e2 )))
             newpos <- unionPositions( positions( e1 ), positions( e2 ) )
	     e1 <- align( e1, newpos, how = "NA" )
	     e2 <- align( e2, newpos, how = "NA" )

	     e1@data <- asSeriesData( callGeneric( e1@data, e2@data ))
	     e1
	   })

setMethod("shift", signature(x="ANY"),
   function(x, k = 1)
      lag(x,  - k)
)

setMethod("shift", "seriesVirtual",
   function(x, k=1)
   {
     if(k != round(k))
     {
       k <- round(k)
       warning("k is not an integer")
     }

     positions(x) <- shiftPositions(positions(x), k)
     x
   })

## a better version of diff for timeSeries is old-style diff.timeSeries() available in FinMetrics

setMethod("diff", "seriesVirtual", function( x, ... ) {
  diff.seriesVirtual <- function(x, lag = 1, differences = 1)
  {
    if(lag < 1 || lag != round(lag))
      stop("lag must be a positive integer")
    if(differences < 1 || differences != round(differences))
      stop("differences must be a positive integer")
  
    newpos <- positions(x)
    nr <- length(newpos)
    nless <- lag * differences
  
    if( nless >= nr )
      return( x[0,] )
  
    newdat <- seriesData(x)
  
    newdat <- diff(newdat, lag, differences)
    # nrow( newpos ) <- nrow( newdat )
    
    x@positions <- newpos[(nr-numRows(newdat)+1):nr]
    x@data <- newdat
    x
  }
  diff.seriesVirtual(x, ...)
})

setMethod("c", signature(x="seriesVirtual"),
   function(x, ...)
   {

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
      newpos <- c(positions(x), positions(y))
      xdat <- seriesData(x)
      ydat <- seriesData(y)
      xvec <- is.null(dim(xdat))
      yvec <- is.null(dim(ydat))
     if( xvec && yvec )
	newdat <- c(xdat, ydat )
     else {
       if( xvec ) 
         xdat <- as.data.frame( xdat )
       if( yvec ) 
         ydat <- as.data.frame( ydat )
       newdat <- rbind( xdat, ydat )
     }
     x@positions <- newpos
     x@data <- newdat
     x
    }
  })

setMethod("cor", signature( x = "seriesVirtual", y = "ANY" ),
  function(x, y = NULL, use = "everything",
         method = c("pearson", "kendall", "spearman"))
  {
    x <- x@data
    callGeneric()
  })

setMethod("cor", signature( x = "seriesVirtual", y = "seriesVirtual" ),
  function(x, y = NULL, use = "everything",
         method = c("pearson", "kendall", "spearman"))
  {
    x <- x@data
    y <- y@data
    callGeneric()
  })

setMethod("cor", signature( x = "ANY", y = "seriesVirtual" ),
  function(x, y = NULL, use = "everything",
         method = c("pearson", "kendall", "spearman"))
  {
    y <- y@data
    callGeneric()
  })


setMethod("var", signature( x = "seriesVirtual" ),
  function(x, y, na.rm = FALSE, use) 
  {
    x <- x@data
    callGeneric()
  })


setMethod( "mean", signature( x = "seriesVirtual" ),
	   function(x, trim = 0.0, na.rm = FALSE, weights = NULL)
	   {
	     x <- x@data
	     callGeneric()
	   })

setMethod( "median", signature( x = "seriesVirtual" ),
	   function(x, na.rm = FALSE)
	   {
	     x <- x@data
	     callGeneric()
	   })

setMethod( "quantile", signature( x = "seriesVirtual" ),
	   function(x, probs = 0:4/4, na.rm = FALSE, ...)
	   {
	     x <- x@data
	     callGeneric()
	   })

setMethod( "aggregate", "seriesVirtual",
           function( x, ... ) aggregateSeries( x, ... ),
            )

setMethod("units", "seriesVirtual", function(x) x@units)

setMethod("deltat", "seriesVirtual",
  function(x, ...)
  {
    if( numRows(x) > 1 ) {
      dfs <- diff( as( x@positions, "numeric"))
      if( all( abs( dfs - dfs[1] ) <= timeDateOptions("sequence.tol")[[1]]))
	return( dfs[1] )
      else
	stop("irregular series -- cannot calculate deltat" )
    }
    return(1)
  })


setMethod("frequency", "seriesVirtual", function(x, ...) 1/deltat(x),
   )

setMethod("window", "seriesVirtual",
  function(x, start=NULL, end=NULL)
  {
    keep <- TRUE
    pos <- positions(x)
    if( !is.null( start ))
      keep <- keep & pos >= start
    if( !is.null( end ))
      keep <- keep & pos <= end
    x[keep,]
  })

##
## make it consistent when you do deltat() on cts objects in SV3
## -- always assuming a period of one year for calendar data as declared in cts help file
## -- other words for a time sequence in one of SV4 time series classes, following holds:
## 	"days"          365       
## 	"weeks"         52 (not enforced)      
##	"months"        12      
##	"quarters"      4      
##	"years"         1 
##
setMethod("deltat", "timeSeries", function(x, ...) {
  mypos <- x@positions
  ## if it's a sequence, the "by" is what we want... 
  ##        if( is( mypos, "timeSequence" ) && !is.missing( mypos@by )) {
  ##            return( mypos@by )
  ##        }
  len <- length(mypos)
  if (numRows(x) <= 1) return(1)
  dfs <- diff(mypos)
  numdfs <- as( dfs, "numeric")
  if( all( abs( numdfs - numdfs[1] ) <= timeDateOptions("sequence.tol")[[1]])) {
    return( dfs[1]/365 )
  } else {
    ## Maybe it's a monthly or yearly sequence
    ## if hours/minutes/seconds all same
    hms1 <- hms(mypos)
    if( any(c( hms1$hour[-1] - hms1$hour[ -len ],
              hms1$minute[-1] - hms1$minute[ -len ],
              hms1$second[-1] - hms1$second[ -len ],
              hms1$ms[-1] - hms1$ms[ -len ] )))			
      stop("irregular series -- cannot calculate deltat" )
    
    ## check if all dates are month end dates
    mdy1 <- mdy(mypos)
    if( !all(is.monthend(mypos )) && any(mdy1$day[-1] - mdy1$day[ -len ]))
      stop( "irregular series -- cannot calculate deltat" )
    mdiffs <- mdy1$month[-1] - mdy1$month[-len] +
      12 * ( mdy1$year[-1] - mdy1$year[-len] )
    if( any( mdiffs - mdiffs[1] ))
      stop( "irregular series -- cannot calculate deltat" )
    ans <- mdiffs[1]/12
    return(ans)
  }
 })

setMethod("sort.list", signature(x="seriesVirtual"),
	function(x, partial=NULL, na.last=TRUE) {
	x <- x@positions
	callGeneric()
	})
	
setMethod("sort", signature(x="seriesVirtual"),
	function(x, partial=NULL, na.last=TRUE) {
		s1 <- sort.list(x, partial=partial, na.last=na.last)
		x[s1,]
	})
	
setMethod( "logb", signature( x = "seriesVirtual" ),
function(x, base = exp(1)){
  x@data = logb(x@data, base=base)
  return(x)
})

setMethod("rep", signature(x = "seriesVirtual"),
function(x, ...)
{
  x@data = rep(x@data, ...)
  x@positions = rep(x@positions, ...)
  x
})
