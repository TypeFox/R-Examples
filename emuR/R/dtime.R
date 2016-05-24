##' time signal times
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export dtime
dtime <- function(dataset, times, single = TRUE, average = TRUE) {
  if(!is.matrix(dataset$data))
    dataset$data <- cbind(dataset$data)
  
  if(!is.matrix(dataset$index)) {
    dataset$index <- rbind(dataset$index)
    dataset$ftime <- rbind(dataset$ftime)
  }
  
  mat <- NULL
  
  for(j in 1:length(times)) {
    left <- dataset$index[j, 1]
    right <- dataset$index[j, 2]
    dat <- dataset$data[left:right,  ]
    if(!is.matrix(dat))
      dat <- cbind(dat)
    lval <- right - left + 1
    left.time <- dataset$ftime[j, 1]
    right.time <- dataset$ftime[j, 2]
    seq.times <- seq(left.time, right.time, length = lval)
    cval <- closest(seq.times, times[j])
    if(single) {
      if(length(cval) > 1) {
        if(average) {
          cval <- mean(cval)
        }
        else {
          cval <- cval[1]
        }
      }
    }
    mat <- rbind(mat, dat[cval,  ])
  }
  
  if(ncol(mat) == 1) {
    c(mat)
  }
  else {
    mat
  }
}














##' Function to extract a vector or matrix from EMU-Trackdata at a single time
##' point of to create another EMU-trackdata object between two times.
##' 
##' A general purpose tool for extracting data from track objects either at a
##' particular time, or between two times. The times can be values in
##' milliseconds or proportional times between zero (the onset) and one (the
##' offset).
##' 
##' This function extracts data from each segment of a trackdata object.
##' 
##' If 'prop=FALSE' the time arguments ('left.time' and 'right.time') are
##' interpreted as millisecond times and each should be a vector with the same
##' length as the number of segments in 'trackdata'.  If 'prop=TRUE' the time
##' arguments should be single values between zero (the onset of the segment)
##' and one (the offset).
##' 
##' If 'right.time' is omitted then a single data point correponding to
##' 'left.time' for each segment is returned.
##' 
##' @aliases dcut dcut.sub
##' @param trackdata An Emu trackdata object.
##' @param left.time Either: a numeric vector of the same length as there are
##' obsverations in trackdata. Or: a single value between 0 and 1. In the first
##' case, the left time boundary of trackdata[n,] is cut at left.time[n], in
##' the second case, and if prop=T, it is cut at that proportional time.
##' @param right.time Either: a numeric vector of the same length as there are
##' obsverations in trackdata. Or: a single value between 0 and 1. In the first
##' case, the right time boundary of trackdata[n,] is cut at right.time[n], in
##' the second case, and if prop=T, it is cut at that proportional time.
##' @param single If TRUE, one value is returned per segment. This applies when
##' the requested time falls between two track frames. When single=TRUE, the
##' preceding value is returned, unless average=TRUE (see below), in which case
##' the average value of the two frames is returned. when the right.time
##' argument is omitted
##' @param average A single element logical vector - see single above. Applies
##' only when the right.times argument is omitted and when single = TRUE
##' @param prop If TRUE left.time and right.time are interpreted as
##' proportions, if FALSE, they are interpreted as millisecond times
##' @return A trackdata object if both 'left.time' and 'right.time' are
##' specified, otherwise a matrix if 'right.time' is unspecified and the
##' trackdata object has multiple columns of data or a vector if right.time' is
##' unspecified and the trackdata object has a single column of data.
##' @author Jonathan Harrington
##' @seealso \code{\link{get_trackdata}}, \code{\link{dplot}}, \code{\link{eplot}}
##' @keywords datagen
##' @examples
##' 
##'      # the data values of the trackdata object at the temporal midpoint
##'      # (midvals is matrix of F1 and F2 data)
##'      dip.fdat[1:10]
##'      midvals <- dcut(dip.fdat, 0.5, prop=TRUE)
##'      midvals[1:10,]
##'      
##'      
##'      # the data values of the trackdata object between 
##'      # extending from 20
##'      # (bet is a trackdata object of F1 and F2 values)
##'      bet <- dcut(dip.fdat, 0.2, 0.8, prop=TRUE)
##'      bet[1]
##'      
##' 
##'      # the data values of the trackdata object at 30 ms after
##'      # the start time of the trackdata object
##'      # (time30 is a matrix of F1 and F2 data
##'      times <- dip.fdat$ftime[,1]+30
##'      times[1:10]
##'      time30 <- dcut(dip.fdat, times)
##'      time30[1:10]
##'      
##' 
##'      # the data values of the trackdata object 
##'      # between the  start time and 30 ms after the start  time
##'      # (int is a trackdata object of F1 and F2 values extending
##'      # from the start of the diphthongs up to 30 ms after the diphthongs)
##'      int <- dcut(dip.fdat, dip.fdat$ftime[,1], times)
##'      int[1]
##'  
##' @export dcut
"dcut" <- function (trackdata, left.time, right.time, 
                    single = TRUE, average = TRUE, prop = FALSE) 
{
  if (prop) {
    if (missing(right.time)) 
      omat <- dextract(trackdata, left.time)
    else omat <- dextract(trackdata, left.time, right.time)
  }
  else {
    if (missing(right.time)) 
      omat <- dtime(trackdata, left.time, single = single, 
                    average = average)
    else {
      if (length(left.time) != nrow(trackdata$ftime)) {
        stop("different number of elements in left.time and $ftime")
      }
      if (length(right.time) != nrow(trackdata$ftime)) {
        stop("different number of elements in right.time and $ftime")
      }
      if (any(left.time < trackdata$ftime[, 1])) 
        stop("some $ftime[,1] values are less than left.time ")
      if (any(right.time > trackdata$ftime[, 2])) 
        stop("some $ftime[,2] values are greater than right.time ")
      if (any(right.time <= left.time)) 
        stop("some right.time values are before the corresponding left.time")
      lval <- nrow(trackdata$index)
      for (j in 1:lval) {
        tdat <- dcut.sub(trackdata[j], left.time[j], 
                         right.time[j])
        if (j == 1) 
          omat <- tdat
        else omat <- bind(omat, tdat)
      }
    }
  }
  if(is.spectral(trackdata$data))
  {
    if(is.trackdata(omat))
    {
      attr(omat$data, "fs") <- attr(trackdata$data, "fs")
      if(!is.spectral(omat$data))
        class(omat$data) <- c(class(omat$data), "spectral")
    }
    else
    {
      attr(omat, "fs") <- attr(trackdata$data, "fs")
      if(!is.spectral(omat))
        class(omat) <- c(class(omat), "spectral")
    }
    
  }
  return(omat)
}



##' @export
"dcut.sub" <- function(trackdata, left.time, right.time)
{
  vals <- trackdata$data
  left <- trackdata$ftime[1]
  right <- trackdata$ftime[2]
  
  if(is.matrix(vals))
    N <- nrow(vals)
  else
    N <- length(vals)
  
  times <- seq(left, right, length = N)
  first <- closest(times, left.time)
  
  if(length(first) > 1)
    first <- round(mean(first))
  
  second <- closest(times, right.time)
  
  if(length(second) > 1)
    second <- round(mean(second))
  
  if(is.matrix(vals))
    trackdata$data <- vals[first:second,  ]
  else
    trackdata$data <- cbind(vals[first:second]
    )
  trackdata$ftime <- cbind(times[first], times[second])
  trackdata$index <- cbind(1, length(first:second))
  
  as.trackdata(trackdata$data, trackdata$index, trackdata$ftime)
}
