##' Calculate gradient summary information for trackdata
##' 
##' Calculates a number of summary measures for a trackdata object: duration,
##' start and end data points, delta values and slope.
##' 
##' \code{track.gradinfo} calculates a number of summary measure for the
##' segments within a trackdata object.  These are useful for data such as
##' kinematic measures where segments might correspond to articulatory
##' movements etc.
##' 
##' Measures returned are: duration, start and end data values (ie. the first
##' and last rows of data for each segment), delta (the difference between the
##' first and last rows of data) and slope (delta divided by the duration).
##' 
##' @param trackdata An Emu trackdata object as returned by
##' \code{\link{get_trackdata}}
##' @return A data frame with one row per segment and columns:
##' \item{duration}{Segment} \item{startN }{The starting value for each segment
##' (start1 is the starting value for the first column) } \item{endN }{The
##' ending value for each segment } \item{deltaN }{The delta value for each
##' segment} \item{slopeN }{The slope value for each segment}
##' 
##' Since the result is a data frame, the columns can be referred to by name
##' (\code{result$duration}) or as matrix columns (\code{result[,1]}).
##' @author Steve Cassidy
##' @seealso \code{\link{get_trackdata}}, \code{\link{dapply}}
##' @keywords misc
##' @examples
##' 
##' data(vowlax)
##' segs = vowlax
##' ## fm has 4 columns
##' data.fm <-vowlax.fdat
##' ## F0 has one
##' data.F0 <- vowlax.fund
##' ## info.fm will have duration, 4xstart, 4xend, 4xdelta, 4xslope
##' info.fm <- track.gradinfo(data.fm)
##' ## this should be true
##' ncol(info.fm) == 1+4+4+4+4
##' 
##' ## info.F0 will have one of each
##' info.F0 <- track.gradinfo(data.F0)
##' ## this should be true
##' ncol(info.F0) == 1+1+1+1+1
##' 
##' ## plot the durations vs delta of the first formant
##' plot(info.F0$duration, info.fm$delta1, type="n", xlab="Duration", ylab="Delta")
##' text(info.fm$duration, info.fm$delta1, labels=label(segs))
##' 
##' ## extract just the delta values from the formant info
##' ## You need to eyeball the data to work out which columns to select
##' delta.fm <- info.fm[,10:13]
##' 
##' @export track.gradinfo
track.gradinfo <- function( trackdata ) {
  ## track.gradinfo  --
  ## generate various bits of information about a trackdata
  ## object:
  ##    - duration
  ##    - start, end: data values at the start and end of the segment
  ##    - delta: the difference between start and end data points
  ##    - slope: the slope of the data (delta/duration)
  ##
  
  result <- dapply(trackdata, track.gradinfo.sub)
  ## all we want is the data which will be one row per segment
  result <- data.frame( result$data )
  
  ## Put appropriate column headers on the data frame
  ## 
  ## this would be better off in track.gradinfo.sub but
  ## because of the way that dapply works it has to go here
  w <- ncol(trackdata$data)
  names(result) <- c("duration", 
                     paste("start", 1:w, sep=""),
                     paste("end", 1:w, sep=""),
                     paste("delta", 1:w, sep=""),
                     paste("slope", 1:w, sep="") )
  return( result )
}









##' track gradinfo sub
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export track.gradinfo.sub
track.gradinfo.sub <- function( data, ftime ) {
  ## track.gradinfo.sub -- 
  ## do the work of track.gradinfo, return the various
  ## measures in the right form for dapply
  n <- nrow(data)
  dur <- ftime[2]-ftime[1]
  ## delta is the difference between the start and end data points
  delta <- diff( data[c(1,n),] )
  ## slope is the delta/duration
  slope <- delta/dur
  data <- matrix( c( dur, data[1,], data[n,], delta, slope ), nrow=1)
  
  ## ftime will be discarded anyway but let's do the right thing
  ## and set the start and end to the segment mid point
  mid <- ftime[1]+dur/2
  ftime <- c(mid, mid)
  return( list( data=data, ftime=ftime ) )
}

