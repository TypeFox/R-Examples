#########################################################
## Methods that define operations on the class "trackdata"
## see also track and frames
#########################################################


##' print trackdata
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export
"print.trackdata"<- function(x, ...)
{
  if(is.null(x$trackname)) 
    cat("trackdata from unknown track.\n")
  else
    cat("trackdata from track:", x$trackname,"\n")
  
  cat("index:\n")
  print(x$index, ...)
  cat("ftime:\n")
  print(x$ftime, ...)
  cat("data:\n")
  print(x$data, ...)
}



##' Expand trackdata
##' 
##' see function
##' 
##' 
##' @aliases [.trackdata
##' @keywords internal
##' @export
"[.trackdata" <- function (dataset, i, j, ...) 
{
  
  
  if (missing(i)) {
    i <- 1:nrow(dataset$index)
  }
  
  
  ftime <- dataset$ftime[i, , drop = FALSE]
  index <- dataset$index[i, , drop = FALSE]
  
  
  datarows <- NULL
  for (ind in 1:nrow(index)) {
    datarows <- c(datarows, seq(from = index[ind, 1], to = index[ind, 
                                                                 2]))
  }
  if (is.matrix(dataset$data)) {
    if (missing(j)) 
      data <- dataset$data[datarows, , drop = FALSE]
    else data <- dataset$data[datarows, j, drop = FALSE]
  } else {
    data <- dataset$data[datarows, drop = FALSE]
  }
  lval <- index[, 2] - index[, 1] + 1
  right <- cumsum(lval)
  left <- right + 1
  left <- left[-length(left)]
  left <- c(1, left)
  nindex <- cbind(left, right)
  dataset$index <- nindex
  dataset$ftime <- ftime
  dataset$data <- data
  return(dataset)
}
















##' summary trackdata
##' 
##' summarizes trackdata objects
##' 
##' 
##' @param object track data object
##' @param \dots see summary
##' @keywords internal
##' @method summary trackdata
##' @export
"summary.trackdata" <- function(object, ...)
{
  if( is.matrix(object$data)){
    dimens <- ncol(object$data)
    len <- nrow(object$data)
  }
  else {
    dimens <- 1
    len <- length(object$data)
  }
  cat("Emu track data from", nrow(object$index), "segments\n\n")
  cat("Data is ", dimens, "dimensional from track", 
      object$trackname,"\n")
  cat("Mean data length is ", len/nrow(object$index), " samples\n")
  invisible()
}










##' Create an Emu trackdata object
##' 
##' Create an Emu trackdata object from a raw data matrix.
##' 
##' Emu trackdata objects contain possibly multi-column numerical data
##' corresponding to a set of segments from a database.  Data for each segment
##' takes up a number of rows in the main \code{data} matrix, the start and end
##' rows are stored in the \code{index} component.  The \code{ftime} component
##' contains the start and end times of the segment data.
##' 
##' Trackdata objects are returned by the \code{\link{get_trackdata}} function.
##' 
##' @param data A two dimensional matrix of numerical data.
##' @param index Segment index, one row per segment, two columns give the start
##' and end rows in the \code{data} matrix for each segment.
##' @param ftime A two column matrix with one row per segment, gives the start
##' and end times in milliseconds for each segment.
##' @param trackname The name of the track.
##' @return The components are bound into a trackdata object.
##' @seealso \code{\link{get_trackdata}} \code{\link{dplot}}
##' @keywords misc
##' @examples
##' 
##' 
##' # make a trackdata object of two data segments
##' data1 <- matrix( 1:10, ncol=2 )
##' data2 <- matrix( 11:20, ncol=2 )
##' 
##' nd1 <- nrow(data1)
##' nd2 <- nrow(data2)
##' index <- rbind( c( 1, nd1 ), c(nd1+1,nd1+nd2) )
##' 
##' times <- rbind( c( 100.0, 110.0 ), c( 200.0, 210.0 ) )
##' 
##' tdata <- as.trackdata( rbind( data1, data2 ), index, times, trackname="fake")
##' 
##' # describe the data
##' summary(tdata)
##' # get the data for the first segment
##' tdata[1]
##' # and the second
##' tdata[2]
##' 
##' 
##' @export as.trackdata
"as.trackdata" <- function( data, index, ftime, trackname="" )
{
  mat <- list( data=as.matrix(data), 
               index=index, 
               ftime=ftime,
               trackname=trackname)
  if( version$major >= 5 ) {
    oldClass(mat) <- "trackdata"
  } else {
    class(mat) <- "trackdata"
  }
  mat
}









##' Test whether an object is an Emu trackdata object
##' 
##' Test whether an object is an Emu trackdata object
##' 
##' 
##' @param object A data object to be tested
##' @return Returns TRUE if the argument is a trackdata object.
##' @seealso \code{\link{get_trackdata}}
##' @keywords misc
##' @export is.trackdata
"is.trackdata" <- function (object) 
{
  return(inherits(object, "trackdata"))
}









##' Produces time-series plots from trackdata
##' 
##' The function produces a plot as a function of time for a single segment or
##' multiple plots as a function of time for several segments.
##' 
##' The function plots a single segment of trackdata as a function of time. If
##' the segment contains multiple tracks, then these will be overlaid. If there
##' are several temporally non-contiguous segments in the trackdata object,
##' each segment is plotted in a different panel by specifying contig=F. This
##' function is not suitable for overlaying trackdata from more than one
##' segments on the same plot as a function of time: for this use dplot().
##' 
##' @param x A trackdata object.
##' @param timestart A single valued numeric vector for setting the time at
##' which the trackdata should start. Defaults to NULL which means that the
##' start time is taken from start(trackdata), i.e. the time at which the
##' trackdata object starts.
##' @param xlim A numeric vector of two values for specifying the time interval
##' over which the trackdata is to be plotted. Defaults to NULL which means
##' that the trackdata object is plotted between between the start time of the
##' first segment and the end time of the last segment.
##' @param ylim Specify a yaxis range.
##' @param labels A character vector the same length as the number of segments
##' in the trackdata object. Each label is plotted at side = 3 on the plotted
##' at the temporal midpoint of each segment in the trackdata object. Defaults
##' to NULL (plot no labels). Labels will only be plotted if xlim=NULL.
##' @param col A single element logical vector. Defaults to T to plot each
##' label type in a different colour
##' @param lty A single element logical vector. Defaults to F.  If TRUE, plot
##' each label type in a different linetype
##' @param type Specify the type of plot. See \link{plot} for the various
##' possibilities
##' @param pch The symbol types to be used for plotting. Should be specified as
##' a numeric vector of the same length as there are unique label classes
##' @param contig A single valued logical vector T or F. If T, then all the
##' segments of the trackdata object are assumed to be temporally contiguous,
##' i.e. the boundaries of the segments are abutting in time and the start time
##' of segment[j,] is the end time of segment[j-1,]. In this case, all the
##' segments of the trackdata object are plotted on the same plot as a function
##' of time. An example of a contiguous trackdata object is coutts.sam. contig
##' = FALSE is when a trackdata object is non-contiguous e.g. all "i:" vowels
##' in a database. An example of a non-contiguous trackdata object is
##' vowlax.fdat. If contig=F then each segment of the trackdata object is
##' plotted separately.
##' @param ...  the same graphical parameters can be supplied to this function
##' as for plot e.g type="l", lty=2 etc.
##' @author Jonathan Harrington
##' @seealso \code{\link{plot}}, \code{\link{dplot}}
##' @keywords dplot
##' @examples
##' 
##' 
##' # a single segment of trackdata (F1) plotted as a function of time.
##' plot(vowlax.fdat[1,1])
##' 
##' # as above, but limits are set for the time axis.
##' plot(vowlax.fdat[1,1], xlim=c(880, 920))
##' 
##' # the the start-time of the x-axis is set to 0 ms, plot F1 and F3, lineplot
##' plot(vowlax.fdat[1,c(1,3)],  timestart=0, type="l")
##' 
##' 
##' # plot F1-F4, same colour, same plotting symbol, between 900 
##' # and 920 ms, type is line and points plot, different linetype per track, no box
##' plot(vowlax.fdat[1,], col="blue", pch=20, xlim=c(900, 920), type="b", lty=TRUE, bty="n")
##' 
##' 
##' # F1 and F2 of six vowels with labels, separate windows
##' par(mfrow=c(2,3))
##' plot(vowlax.fdat[1:6,1:2], contig=FALSE, labels=vowlax.l[1:6], ylab="F1 and F2", 
##' xlab="Time (ms)", type="b", ylim=c(300, 2400))
##' 
##' # As above, timestart set to zero, colour set to blue, different plotting
##' # symbols for the two tracks
##' plot(vowlax.fdat[1:6,1:2], contig=FALSE, labels=vowlax.l[1:6], ylab="F1 and F2", 
##' xlab="Time (ms)", type="b", col="blue", pch=c(1,2),  ylim=c(300, 2400), timestart=0)
##' 
##' # RMS energy for the utterance 'just relax said Coutts'
##'  plot(coutts.rms, type="l")
##' # as above a different colour
##'  plot(coutts.rms, type="l", col="pink")
##' # as above, linetype 2, double line thickness, no box, times reset to 0 ms
##'  plot(coutts.rms, type="l", col="pink", lty=2, lwd=2, bty="n", timestart=0)
##' # as above but plotted as non-contiguous segments, i.e one segment per panel
##'  par(mfrow=c(2,3))
##'  plot(coutts.rms, type="l", col="pink", lty=2, lwd=2, bty="n", timestart=0, contig=FALSE)
##' # plot with labels
##'  labels = label(coutts)
##' par(mfrow=c(1,1))
##'  plot(coutts.rms, labels=labels, type="l", bty="n")
##' # as above, double line-thickness, green, line type 3, no box, 
##' # time start 0 ms with x and y axis labels
##'  plot(coutts.rms, labels=labels, type="l", lwd=2, 
##'       col="green", lty=3,  bty="n", timestart=0, xlab="Time (ms)", ylab="Amplitude")
##' # as above with a different plotting symbol for the points
##' par(mfrow=c(2,3))
##'  plot(coutts.rms, labels=labels, type="b", lwd=2, col="green", 
##'       timestart=0, bty="n", contig=FALSE, pch=20)
##'  
##' 
##' 
##' @export
`plot.trackdata` <- function (x, timestart = NULL, xlim = NULL, 
                              ylim = NULL, labels = NULL, col = TRUE, 
                              lty = FALSE, type="p", pch=NULL, 
                              contig = TRUE, ...) 
{
  trackdata <- x
  N <- nrow(trackdata$data)
  if(is.logical(col))
  {
    if (col) 
      col <- 1:ncol(trackdata)
    else
      col <- rep(1, ncol(trackdata))
  }
  else
  {
    if(length(col)!=ncol(trackdata))
      col <- rep(col[1], ncol(trackdata))
  }
  
  if(is.logical(lty))
  {
    if (lty) 
      lty <- 1:ncol(trackdata)
    else
      lty <- rep(1, ncol(trackdata))
  }
  else
  {
    if(length(lty)!=ncol(trackdata))
      lty <- rep(lty[1], ncol(trackdata))
  }
  if(is.null(pch))
    pch <- rep(1, ncol(trackdata))
  else
  {
    if(length(pch)!=ncol(trackdata))
      pch <- rep(pch[1], ncol(trackdata))
  }
  
  
  n <- nrow(trackdata)
  if (!is.null(xlim)) 
    labels <- NULL
  if (!is.null(labels)) {
    if (length(labels) != nrow(trackdata)) 
      stop("if labels are supplied, there must be one label per segment")
    label.times <- apply(trackdata$ftime, 1, mean)
    boundary.times <- c(trackdata$ftime[, 1], trackdata$ftime[n])
  }
  if (n > 1 & contig) {
    inds <- cbind(1, N)
    ftime <- cbind(trackdata$ftime[1, 1], trackdata$ftime[n, 
                                                          2])
    trackdata <- as.trackdata(trackdata$data, inds, ftime)
  }
  if (!is.null(xlim)) {
    if (nrow(trackdata) != 1) 
      stop("can't specify xlim if there's more than one segment")
  }
  left <- trackdata$ftime[1]
  right <- trackdata$ftime[2]
  times <- seq(left, right, length = nrow(trackdata$data))
  if (!is.null(timestart)) {
    times <- times - left + timestart
    if (!is.null(labels)) {
      label.times <- label.times - left + timestart
      boundary.times <- boundary.times - left + timestart
    }
  }
  data <- trackdata$data
  if (nrow(trackdata) == 1) {
    if (is.null(xlim)) 
      xlim <- range(times)
    if (is.null(ylim)) 
      ylim <- range(data)
    for (k in 1:ncol(data)) {
      if(k==ncol(data))
        graphics::plot(times, data[, k], xlim = xlim, ylim = ylim, 
             col = col[k], lty = lty[k], pch=pch[k], type=type, ...)
      else
        graphics::plot(times, data[, k], xlim = xlim, ylim = ylim, 
             col = col[k], lty = lty[k], pch=pch[k], xlab="", ylab="", main="", axes=FALSE, bty="n", type=type)
      graphics::par(new = TRUE)
    }
    graphics::par(new = FALSE)
    if (!is.null(labels)) {
      if (length(boundary.times) > 2) 
        graphics::abline(v = boundary.times)
      graphics::mtext(labels, at = label.times)
    }
  }
  else {
    if (is.null(labels)) 
      labels <- rep("", nrow(trackdata))
    for (j in 1:nrow(trackdata)) {
      graphics::plot(trackdata[j, ], timestart = timestart, xlim = xlim, 
           ylim = ylim, labels = labels[j], col=col, lty=lty, type=type, pch=pch,contig = TRUE, ...)
    }
  }
}



##' @export
"bark.trackdata" <- function(f, ...)
{
  trackdata = f
  if(is.spectral(trackdata$data))
    return(bark.spectral(trackdata))
  else
  {
    trackdata$data <- bark(trackdata$data)
    return(trackdata)
  }
}


##' @export
"mel.trackdata" <- function(a)
{
  trackdata = a
  if(is.spectral(trackdata$data))
    return(mel.spectral(trackdata))
  else
  {
    trackdata$data <- mel(trackdata$data)
    return(trackdata)
  }
}













##' function to find the frequencies of a spectral object
##' 
##' Find the frequencies of a spectral object.
##' 
##' 
##' @param specdata A spectral object
##' @return A vector of the frequencies at which the columns of a spectral
##' matrix occur.
##' @author Jonathan Harrington
##' @keywords attribute
##' @examples
##' 
##' trackfreq(vowlax.dft.5)
##' # Frequency components between 1000 and 2000 Hz
##' trackfreq(vowlax.dft.5[,1000:2000])
##' # All frequency components of a trackdata object except the d.c. offset
##' trackfreq(fric.dft[,-1])
##' # All frequency components except the d.c. offset
##' # and except frequencies above 5000 Hz
##' trackfreq(fric.dft[,-c(1, 5000:20000)])
##' # Note the following syntax if the spectral object is a vector
##' # Frequencies 1000-3000 Hz
##' trackfreq(e.dft[1000:3000])
##' 
##' 
##' 
##' @export trackfreq
"trackfreq" <- function(specdata){
  if(is.trackdata(specdata))
    return(attr(specdata$data, "fs"))
  else
    return(attr(specdata, "fs"))
}










##' get trackkeywrd
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export get.trackkeywrd
"get.trackkeywrd" <- function (fname) 
{
  line <- readLines(fname, n = 2)
  if (length(line) < 2) {
    return(NULL)
  }
  
  line <- splitstring(line[2], " ")
  if ((length(line) == 3) && (line[2] == "Trackname")) {
    trackname <- line[3]
  }
  else {
    return(NULL)
  }
  if (trackname  !=  "") {
    return(trackname)
  }
  else {
    return(NULL)
  }
}










##' Duration of trackdata elements
##' 
##' Duration of segments is calculated for each element in the trackdata object
##' 
##' 
##' @param x a trackdata object
##' @return a vector of durations
##' @author Jonathan Harrington
##' @keywords internal
##' @export
"dur.trackdata" <- function (x) 
{
  x$ftime[,2] - x$ftime[,1]
}









##' frames
##' 
##' Get frames from trackdata objects
##' 
##' 
##' @param trackdata an object of class trackdata
##' @return Data frames from the input object.
##' @author Jonathan Harrington
##' @seealso \code{\link{trackdata}}
##' @keywords utilities
##' @export frames
"frames" <- function(trackdata)
{
  if(!(is.trackdata(trackdata)))
    stop ("Object must be of class trackdata")
  trackdata$data
}
