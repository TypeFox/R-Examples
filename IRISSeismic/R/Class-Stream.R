##
##    S4 classes for handling lists of seismic traces
##
##    Copyright (C) 2012  Mazama Science, Inc.
##    by Jonathan Callahan, jonathan@mazamascience.com
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


################################################################################
# R class for a Stream obj.
#
# This is a port of some of the functionality found in obspy.core.stream.py 
# 
#   http://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html
#
################################################################################

################################################################################
# Class Stream
#
# A Stream contains a list of Trace objects.
# The 'url' slot contains the URL of the webservice request that returned the 
# stream -- useful for debugging.
#
################################################################################

setClass("Stream", 
  # typed slots (aka attributes) for class Stream
  representation(url = "character",
                 requestedStarttime = "POSIXct",
                 requestedEndtime = "POSIXct",
                 act_flags = "integer",
                 io_flags = "integer",
                 dq_flags = "integer",
                 timing_qual = "numeric",
                 traces = "list"),
  # default values for slots
  prototype(url = "",
            requestedStarttime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
            requestedEndtime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
            act_flags = rep(as.integer(0),8),
            io_flags = rep(as.integer(0),8),
            dq_flags = rep(as.integer(0),8),
            timing_qual = as.numeric(NA),
            traces = list(new("Trace")))
)

################################################################################
# Basic methods to operate directly on a Stream object's Trace data.
#
# The 'sd' method returns a value for each Trace in the Stream.
#
# In the obspy version of the code, the standard deviation method is called 'std'.
# In this package, the name is changed to 'sd' to match the name for the 
# existing R function for standard deviation.
#
# The 'pmax', 'pmin', etc. methods return a value fore ach Trace in the Stream.
# The 'max', 'mean', etc. methods return a single value for all Traces in the Stream.
################################################################################

# NOTE:  The 'lapply' function returns the list that results from applying the
# NOTE:  in-line function to each element of x@traces.
# NOTE:  The 'unlist' function converts the list returned by 'lapply' into 
# NOTE:  a numeric vector.

# Return unique Ids ------------------------------------------------------------

if (!isGeneric("uniqueIds")) {
  setGeneric("uniqueIds", function(x) {
    standardGeneric("uniqueIds")
  })
}
uniqueIds.Stream <- function(x, na.rm=FALSE) {
  ids <- unlist(lapply(x@traces, slot, "id"))
  return( unique(ids) )
}
setMethod("uniqueIds", signature(x="Stream"), function(x) uniqueIds.Stream(x))


# Parallel Length --------------------------------------------------------------

if (!isGeneric("parallelLength")) {
  setGeneric("parallelLength", function(x) {
    standardGeneric("parallelLength")
  })
}
parallelLength.Stream <- function(x, na.rm=FALSE) {
  return( unlist(lapply(x@traces, function(element) length(element))) )
}
setMethod("parallelLength", signature(x="Stream"), function(x) parallelLength.Stream(x))


# Parallel Maximum -------------------------------------------------------------

if (!isGeneric("parallelMax")) {
  setGeneric("parallelMax", function(x, na.rm) {
    standardGeneric("parallelMax")
  })
}
parallelMax.Stream <- function(x, na.rm) {
  return( unlist(lapply(x@traces, function(element) max(element, na.rm=na.rm))) )
}
setMethod("parallelMax", signature(x="Stream", na.rm="logical"), function(x, na.rm) parallelMax.Stream(x, na.rm=na.rm))
setMethod("parallelMax", signature(x="Stream", na.rm="missing"), function(x, na.rm) parallelMax.Stream(x, na.rm=FALSE))


# Parallel Mean ----------------------------------------------------------------

if (!isGeneric("parallelMean")) {
  setGeneric("parallelMean", function(x, na.rm) {
    standardGeneric("parallelMean")
  })
}
parallelMean.Stream <- function(x, na.rm) {
  return( unlist(lapply(x@traces, function(element) mean(element, na.rm=na.rm))) )
}
setMethod("parallelMean", signature(x="Stream", na.rm="logical"), function(x, na.rm) parallelMean.Stream(x, na.rm=na.rm))
setMethod("parallelMean", signature(x="Stream", na.rm="missing"), function(x, na.rm) parallelMean.Stream(x, na.rm=FALSE))


# Parallel Minimum--------------------------------------------------------------

if (!isGeneric("parallelMin")) {
  setGeneric("parallelMin", function(x, na.rm) {
    standardGeneric("parallelMin")
  })
}
parallelMin.Stream <- function(x, na.rm) {
  return( unlist(lapply(x@traces, function(element) min(element, na.rm=na.rm))) )
}
setMethod("parallelMin", signature(x="Stream", na.rm="logical"), function(x, na.rm) parallelMin.Stream(x, na.rm=na.rm))
setMethod("parallelMin", signature(x="Stream", na.rm="missing"), function(x, na.rm) parallelMin.Stream(x, na.rm=FALSE))


# Parallel Median --------------------------------------------------------------

if (!isGeneric("parallelMedian")) {
  setGeneric("parallelMedian", function(x, na.rm) {
    standardGeneric("parallelMedian")
  })
}
parallelMedian.Stream <- function(x, na.rm) {
  return( unlist(lapply(x@traces, function(element) median(element, na.rm=na.rm))) )
}
setMethod("parallelMedian", signature(x="Stream", na.rm="logical"), function(x, na.rm) parallelMedian.Stream(x, na.rm=na.rm))
setMethod("parallelMedian", signature(x="Stream", na.rm="missing"), function(x, na.rm) parallelMedian.Stream(x, na.rm=FALSE))


# Parallel Standard Deviation --------------------------------------------------

if (!isGeneric("parallelSd")) {
  setGeneric("parallelSd", function(x, na.rm) {
    standardGeneric("parallelSd")
  })
}
parallelSd.Stream <- function(x, na.rm) {
  return( unlist(lapply(x@traces, function(element) sd(element, na.rm=na.rm))) )
}
setMethod("parallelSd", signature(x="Stream", na.rm="logical"), function(x, na.rm) parallelSd.Stream(x, na.rm=na.rm))
setMethod("parallelSd", signature(x="Stream", na.rm="missing"), function(x, na.rm) parallelSd.Stream(x, na.rm=FALSE))


# Parallel RMS -----------------------------------------------------------------

if (!isGeneric("parallelRms")) {
  setGeneric("parallelRms", function(x, na.rm) {
    standardGeneric("parallelRms")
  })
}
parallelRms.Stream <- function(x, na.rm) {
  return( unlist(lapply(x@traces, function(element) rms(element, na.rm=na.rm))) )
}
setMethod("parallelRms", signature(x="Stream", na.rm="logical"), function(x, na.rm) parallelRms.Stream(x, na.rm=na.rm))
setMethod("parallelRms", signature(x="Stream", na.rm="missing"), function(x, na.rm) parallelRms.Stream(x, na.rm=FALSE))


# Parallel RMS variance --------------------------------------------------------

if (!isGeneric("parallelRmsVariance")) {
  setGeneric("parallelRmsVariance", function(x, na.rm) {
    standardGeneric("parallelRmsVariance")
  })
}
parallelRmsVariance.Stream <- function(x, na.rm) {
  return( unlist(lapply(x@traces, function(element) rmsVariance(element, na.rm=na.rm))) )
}
setMethod("parallelRmsVariance", signature(x="Stream", na.rm="logical"), function(x, na.rm) parallelRmsVariance.Stream(x, na.rm=na.rm))
setMethod("parallelRmsVariance", signature(x="Stream", na.rm="missing"), function(x, na.rm) parallelRmsVariance.Stream(x, na.rm=FALSE))


# Length -----------------------------------------------------------------------

length.Stream <- function(x) {
  return( sum(parallelLength(x)) )
}
setMethod("length", signature(x="Stream"), function(x) length.Stream(x))


# Global Maximum ---------------------------------------------------------------

max.Stream <- function(x, ..., na.rm=FALSE) {
  return( max(parallelMax(x, na.rm=na.rm)) )
}
setMethod("max", signature(x="Stream"), function(x, ...) max.Stream(x, ...))


# Global Mean ------------------------------------------------------------------

mean.Stream <- function(x, ...) {
  data <- unlist(lapply(x@traces, slot, "data"))
  return( mean(data, ...) )
}
setMethod("mean", signature(x="Stream"), function(x, ...) mean.Stream(x, ...))


# Global Minimum ---------------------------------------------------------------

min.Stream <- function(x, ..., na.rm=FALSE) {
  return( min(parallelMin(x, na.rm=na.rm)) )
}
setMethod("min", signature(x="Stream"), function(x, ...) min.Stream(x, ...))


# Global Median ----------------------------------------------------------------

median.Stream <- function(x, na.rm) {
  data <- unlist(lapply(x@traces, slot, "data"))
  return( median(data, na.rm=na.rm) )
}
# NOTE:  method signature must match generic signature for function 'median' with arguments: 'x', 'na.rm'
setMethod("median", signature(x="Stream", na.rm="logical"), function(x, na.rm) median.Stream(x, na.rm=na.rm))
setMethod("median", signature(x="Stream", na.rm="missing"), function(x, na.rm) median.Stream(x, na.rm=FALSE))


# Global Standard Deviation ----------------------------------------------------

sd.Stream <- function(x, na.rm) {
  data <- unlist(lapply(x@traces, slot, "data"))
  return( sd(data, na.rm=na.rm) )
}
# NOTE:  method signature must match generic signature for function 'sd' with arguments: 'x', 'na.rm'
setMethod("sd", signature(x="Stream", na.rm="logical"), function(x, na.rm) sd.Stream(x, na.rm=na.rm))
setMethod("sd", signature(x="Stream", na.rm="missing"), function(x, na.rm) sd.Stream(x, na.rm=FALSE))


# Global Root Mean Square ------------------------------------------------------

if (!isGeneric("rms")) {
  setGeneric("rms", function(x, na.rm) {
    standardGeneric("rms")
  })
} 
rms.Stream <- function(x, na.rm) {
  data <- unlist(lapply(x@traces, slot, "data"))
  return( sqrt( mean((data)^2, na.rm=na.rm) ) )
}
setMethod("rms", signature("Stream", na.rm="logical"), function(x, na.rm) rms.Stream(x, na.rm=na.rm))
setMethod("rms", signature("Stream", na.rm="missing"), function(x, na.rm) rms.Stream(x, na.rm=FALSE))


# Global Root Mean Square Variance ---------------------------------------------

if (!isGeneric("rmsVariance")) {
  setGeneric("rmsVariance", function(x, na.rm) {
    standardGeneric("rmsVariance")
  })
} 
rmsVariance.Stream <- function(x, na.rm) {
  data <- unlist(lapply(x@traces, slot, "data"))
  mean <- mean(data, na.rm=na.rm)
  n <- length(data)
  return( sqrt( sum( (data-mean)^2 ) / n ) )
}
setMethod("rmsVariance", signature("Stream", na.rm="logical"), function(x, na.rm) rmsVariance.Stream(x, na.rm=na.rm))
setMethod("rmsVariance", signature("Stream", na.rm="missing"), function(x, na.rm) rmsVariance.Stream(x, na.rm=FALSE))


################################################################################
# Method to return a new Stream object where every Trace has been multiplied
# by a constant.
################################################################################

if (!isGeneric("multiplyBy")) {
  setGeneric("multiplyBy", function(x, y) {
    standardGeneric("multiplyBy")
  })
}
multiplyBy.Stream <- function(x, y) {
  traces <- lapply(x@traces, function(element) multiplyBy(element, y=y))
  return( new("Stream", url=x@url, requestedStarttime=x@requestedStarttime, requestedEndtime=x@requestedEndtime,
              act_flags=x@act_flags, io_flags=x@io_flags, dq_flags=x@dq_flags, timing_qual=x@timing_qual,
              traces=traces) )
}
setMethod("multiplyBy", signature("Stream", y="numeric"), function(x, y) multiplyBy.Stream(x, y=y))


################################################################################
# Method to return a list of all gaps/overlaps of Traces in the Stream object
#
# getGaps(Stream, min_gap)
#   min_gap -- minimum gap (sec) below which gaps will be ignored (default=NULL)
#
################################################################################

if (!isGeneric("getGaps")) {
  setGeneric("getGaps", function(x, min_gap) {
    standardGeneric("getGaps")
  })
}

getGaps.Stream <- function(x, min_gap) {

  # Sanity check -- single SNCL
  num_ids <- length(uniqueIds(x))
  if (num_ids > 1) {
    stop(paste("getGaps.Stream:",num_ids,"unique ids encountered in Stream."))
  }
  
  # NOTE:  The number of potential gaps includes an initial gap, a final gap
  # NOTE:  and a gap between each trace. Thus num_gaps = 2 + num_headers - 1.
  # NOTE:  We will set up an array for num_header+1 gaps but insert 0 if the
  # NOTE:  gap is too small to count.

  # NOTE:  Use difftime() instead of just subtracting to guarantee that units are "secs"
  
  # Extract a list of TraceHeaders from the list of Traces
  headers <- lapply(x@traces, slot, "stats")
  num_headers <- length(headers)
  
  # Sanity check -- single sampling_rate
  sampling_rates <- sapply(headers, slot, "sampling_rate")
  num_rates <- length(unique(round(sampling_rates,digits=4)))
  if (num_rates > 1) {
    stop(paste("getGaps.Stream:",num_rates,"unique sampling rates encountered in Stream."))    
  }
  sampling_rate <- sampling_rates[1]
  
  # Set up arrays for information about gaps/overlaps
  gaps <- numeric(num_headers+1)
  nsamples <- integer(num_headers+1)   
  
  # Set min_gap and make sure it is at least 1/sampling_rate
  min_gap <- ifelse(is.null(min_gap), 1/sampling_rate, min_gap)
  min_gap <- max(min_gap, 1/sampling_rate)

  # NOTE:  The delta we calculate here has a valid datapoint at it's start and end.
  # NOTE:  For the purposes of calculating the correct number of missing points we
  # NOTE:  will calculate delta as the total time between points minus 1/sampling rate.
  # NOTE:  That is how many extra points could be shoved into this gap.
  
  # Initial gap (no overlap possible)
  delta <- as.numeric(difftime(headers[[1]]@starttime, x@requestedStarttime, units="secs")) - 1/sampling_rate
  if (abs(delta) > min_gap) {
    gaps[1] <- delta
    nsamples[1] <- as.integer(round(abs(delta) * sampling_rate))     
  } else {
    gaps[1] <- 0
    nsamples[1] <- 0
  }
  
  # Inter-trace gaps and overlaps
  if (num_headers > 1) {
    for ( i in seq(from=2, to=num_headers) ) {
      h1 <- headers[[i-1]]
      h2 <- headers[[i]]
      delta <- difftime(h2@starttime, h1@endtime, units="secs") - 1/sampling_rate
      if (abs(delta) > min_gap) {
        gaps[i] <- delta
        nsamples[i] <- as.integer(round(abs(delta) * sampling_rate))     
      } else {
        gaps[i] <- 0
        nsamples[i] <- 0
      }
    }    
  }
  
  # Final gap (no overlap possible)
  delta <- as.numeric(difftime(x@requestedEndtime, headers[[num_headers]]@endtime, units="secs")) - 1/sampling_rate
  if (abs(delta) > min_gap) {
    gaps[num_headers+1] <- delta
    nsamples[num_headers+1] <- as.integer(round(abs(delta) * sampling_rate))     
  } else {
    gaps[num_headers+1] <- 0
    nsamples[num_headers+1] <- 0
  }
  
  gap_list <- list(gaps=gaps,
                   nsamples=nsamples)    
    
  return(gap_list)
  
}

# All parameters specified
setMethod("getGaps", signature(x="Stream", min_gap="numeric"), 
          function(x, min_gap) getGaps.Stream(x, min_gap))
# min_gap missing
setMethod("getGaps", signature(x="Stream", min_gap="missing"), 
          function(x, min_gap) getGaps.Stream(x, NULL))


################################################################################
# Method to return a vector of up/down times identifying the start- and end-times
# of Traces in the Stream object.  Care is taken to identify overlapping Traces
# gaps below a minimum threshold and remove those timepoints.
#
# getUpDownTimes(Stream, min_signal, min_gap)
#   min_signal -- minimum signal (sec) below which Traces will be ignored (default=NULL)
#   min_gap -- minimum gap (sec) below which gaps will be ignored (default=NULL)
#
################################################################################

if (!isGeneric("getUpDownTimes")) {
  setGeneric("getUpDownTimes", function(x, min_signal, min_gap) {
    standardGeneric("getUpDownTimes")
  })
}
getUpDownTimes.Stream <- function(x,
                                  min_signal,
                                  min_gap) {
  
  # Extract a list of ids from the list of Traces
  num_ids <- length(uniqueIds(x))
  if (num_ids > 1) {
    stop(paste("getUpDownTimes.Stream:",num_ids,"unique ids encountered in Stream."))
  }
  
  # Extract a list of TraceHeaders
  headerList <- lapply(x@traces, slot, "stats")
  
  # Extract a list of start- and end-times from the list of Traces
  starttimeList <- lapply(headerList, slot, "starttime")
  endtimeList <- lapply(headerList, slot, "endtime")

  # NOTE:  unlist() will convert 'POSIXct' to 'numeric' so we need to go through 'character'
  # NOTE:  This seems inefficient but we're only dealing with a small number of items.
  starttimes <- as.POSIXct(unlist(lapply(starttimeList, strftime, format="%Y-%m-%dT%H:%M:%OS", tz="GMT")), 
                         format="%Y-%m-%dT%H:%M:%OS", tz="GMT")
  endtimes <- as.POSIXct(unlist(lapply(endtimeList, strftime, format="%Y-%m-%dT%H:%M:%OS", tz="GMT")), 
                         format="%Y-%m-%dT%H:%M:%OS", tz="GMT")

  # Determine signal duration for each trace
  signal_durations <- difftime(endtimes, starttimes, units="sec")

  good_traces_flag <- signal_durations >= min_signal 

  # Cull to remove those Traces whose duration is < min_signal  
  headerList <- headerList[good_traces_flag]
  starttimes <- starttimes[good_traces_flag]
  endtimes <- endtimes[good_traces_flag]
  num_headers <- length(headerList)  
  
  # Return if there is only one trace
  if (num_headers == 1) {
    up_down_times <- c(starttimes[1], endtimes[1])    
    return(up_down_times)
  }
  
  # Create a dummy vector of POSIXct values of length 2 X num_headers.
  up_down_times <- c(starttimes,endtimes)
  
  # Now go through and insert appropriate times, checking for any overlaps or spaces < min_gap.
  # NOTE:  up_down_times[1] and up_down_times[num_headers*2] are already in place.
  for (i in seq(num_headers-1)) {
	# get the difference in seconds between the current end time and the next start time.
    delta <- difftime(starttimes[i+1], endtimes[i], units="secs")
    if (delta < 0 || delta < min_gap) {
	  # insert NA empty space to indicate continuity
      up_down_times[(2*i)] <- NA
      up_down_times[(2*i)+1] <- NA
    } else {
	  # insert our current end time and write the next start time to the vector, indicating discontinuity
      up_down_times[(2*i)] <- endtimes[i]
      up_down_times[(2*i)+1] <- starttimes[i+1]
    }
  }
  
  # remove all of the NA values that represent continuity
  # return the vector of start,end,start,end....
  return(stats::na.omit(up_down_times))
    
}
setMethod("getUpDownTimes", signature(x="Stream", min_signal="numeric", min_gap="numeric"), 
          function(x, min_signal, min_gap) getUpDownTimes.Stream(x, min_signal, min_gap))
setMethod("getUpDownTimes", signature(x="Stream", min_signal="missing", min_gap="missing"), 
          function(x, min_signal, min_gap) getUpDownTimes.Stream(x, min_signal=30, min_gap=60))


################################################################################
# Method to return a new Stream object sliced out of an existing Stream
#
# Unlike the ObsPy method, this method does no padding and does not retain empty
# Traces.  The returned Stream will always be a subset of the original Stream
#
# Slicing is rounded to the nearest second.
#
# http://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.slice.html#obspy.core.stream.Stream.slice
################################################################################

if (!isGeneric("slice")) {
  setGeneric("slice", function(x, starttime, endtime) {
    standardGeneric("slice")
  })
} 
slice.Stream <- function(x, starttime, endtime) {
  num_traces <- length(x@traces)
  stream_start <- x@traces[[1]]@stats@starttime
  stream_end <- x@traces[[num_traces]]@stats@endtime
  
  # Sanity check
  if (starttime >= endtime) {
    stop(paste("slice.Stream: requested starttime \"", starttime, "\" >= requested endtime \"", endtime, "\""))    
  }
  if (starttime >= stream_end) {
    stop(paste("slice.Stream: requested starttime \"", starttime, "\" >= Stream endtime \"", stream_end, "\""))
  }
  if (endtime <= stream_start) {
    stop(paste("slice.Stream: requested endtime \"", endtime, "\" <= Stream starttime \"", stream_start, "\""))
  }
  
  # Set up a new list for traces we are keeping
  traces <- list()
  
  # TODO:  Get rid of append() in slice.Stream().
  
  # Check each trace to see if it should be ignored, included entirely or included after slicing
  for (i in seq(num_traces)) {
    tr <- x@traces[[i]]
    if (starttime >= tr@stats@endtime || endtime <= tr@stats@starttime) {
      # skip this trace
    } else if (starttime <= tr@stats@starttime && endtime >= tr@stats@endtime) {
      # this traces does not need to be sliced
      traces <- append(traces,tr)
    } else {
      sliced_trace <- slice(tr, starttime, endtime)
      traces <- append(traces, sliced_trace)
    }
  }
  
  return( new("Stream", url=x@url, requestedStarttime=starttime, requestedEndtime=endtime,
              act_flags=x@act_flags, io_flags=x@io_flags, dq_flags=x@dq_flags, timing_qual=x@timing_qual,
              traces=traces) )
}
setMethod("slice", signature(x="Stream", starttime="POSIXct", endtime="POSIXct"),
          function(x, starttime, endtime) slice.Stream(x, starttime=starttime, endtime=endtime))


################################################################################
# Method to merge Traces in a Stream into a single Trace, replacing gaps with 
# values determined by the fillMethod.
################################################################################

if (!isGeneric("mergeTraces")) {
  setGeneric("mergeTraces", function(x, fillMethod) {
    standardGeneric("mergeTraces")
  })
} 

mergeTraces.Stream <- function(x, fillMethod) {
  
  # Return immediately if is only one trace with no initial or final gap
  if (sum(getGaps(x)$nsamples) == 0) { return(x) }
  
  num_traces <- length(x@traces)
  
  gapInfo <- getGaps(x)
  num_gaps <- length(gapInfo$nsamples)
  
  # Sanity check
  if (num_gaps != num_traces+1) {
    stop(paste("mergeTraces.Stream: num_gaps (",num_gaps,") should be one more than num_traces (",num_traces,")", sep=""))
  }
  
  # NOTE:  Increment totalPoints by one to account for point at the very end.
  # NOTE:  For example, the sequence *--*--*--* has 4 points but only three time intervals.
  totalSecs <- as.numeric(difftime(x@requestedEndtime, x@requestedStarttime, units="secs"))
  totalPoints <- as.integer(round(totalSecs) * x@traces[[1]]@stats@sampling_rate) + 1

  # NOTE:  Setting up a list of vectors that we will concatenate in one fell swoop at the end.
  # NOTE:  This avoids the incremental growth and copying that is incredibly slow.
  # NOTE:  See R Inferno Chapter 2 -- Growing Objects
  # NOTE:    http://www.burns-stat.com/pages/Tutor/R_inferno.pdf
  
  num_vectors <- num_gaps + num_traces
  dataList <- vector('list',num_vectors)
  
  # Fill in the data for gaps and traces
  if (fillMethod == "fillNA") {       
    for (i in seq(num_traces)) {
      dataList[[2*i-1]] <- rep(NA,gapInfo$nsamples[i])
      dataList[[2*i]] <- x@traces[[i]]@data    
    }
    dataList[[num_vectors]] <- rep(NA,gapInfo$nsamples[[num_gaps]])    
  } else if (fillMethod == "fillZero") {    
    for (i in seq(num_traces)) {
      dataList[[2*i-1]] <- rep(0,gapInfo$nsamples[i])
      dataList[[2*i]] <- x@traces[[i]]@data    
    }
    dataList[[num_vectors]] <- rep(0,gapInfo$nsamples[[num_gaps]])    
  } else {
    stop(paste("mergeTraces.Stream: unknown fillMethod '", fillMethod, "'",sep=""))      
  }
  
  # Create a single data vector
  data <- unlist(dataList)
  
  missing_points <- totalPoints - length(data)
  
  # Sanity check -- we hope to be within twice the sampling rate of a complete accounting
  if ( missing_points > ceiling(2 * x@traces[[1]]@stats@sampling_rate) )  {
    stop(paste("mergeTraces.Stream:", missing_points, "unaccounted for points after merge"))
  } else if ( missing_points < ceiling(-2 * x@traces[[1]]@stats@sampling_rate) ) {
    stop(paste("mergeTraces.Stream:", abs(missing_points), "extra points after merge"))
  }
  
  # Guarantee that we return the correct number of points
  if (missing_points > 0) {
    data <- c(data,rep(NA,missing_points))
  }
  
  
  # Create a new TraceHeader
  stats <- x@traces[[1]]@stats
  stats@npts <- as.integer(totalPoints)
  stats@starttime <- x@requestedStarttime
  stats@endtime <- x@requestedEndtime
  stats@processing <- append(stats@processing,paste(num_traces," traces merged into a single trace using method '",fillMethod,"'",sep=""))
  
  # Other Trace info
  id <- x@traces[[1]]@id
  Sensor <- x@traces[[1]]@Sensor
  InstrumentSensitivity <- x@traces[[1]]@InstrumentSensitivity
  InputUnits <- x@traces[[1]]@InputUnits
  
  traces <- list( new("Trace", id, stats, Sensor, InstrumentSensitivity, InputUnits, data=data[1:totalPoints]) )
     
  return( new("Stream", url=x@url, requestedStarttime=x@requestedStarttime, requestedEndtime=x@requestedEndtime,
              act_flags=x@act_flags, io_flags=x@io_flags, dq_flags=x@dq_flags, timing_qual=x@timing_qual,
              traces=traces) )
}

# All parameters specified
setMethod("mergeTraces", signature(x="Stream", fillMethod="character"),
          function(x, fillMethod) mergeTraces.Stream(x, fillMethod=fillMethod))
# fillMethod missing
setMethod("mergeTraces", signature(x="Stream", fillMethod="missing"),
          function(x, fillMethod) mergeTraces.Stream(x, fillMethod="fillNA"))

################################################################################
# Various methods and functions for upDownTimes
#
# upDownTimes is a vector of class POSIXct that contains the GMT times at which
# a channel or channels is on (up) or off (down).  The first and all odd elements
# in the list are datetimes at which data collection starts.  Even numbered
# elements correspond to datetimes at which data collection stops.
################################################################################

################################################################################
# Plotting upDownTimes from Streams or pre-generated vectors of upDownTimes
################################################################################

if (!isGeneric("plotUpDownTimes")) {
  setGeneric("plotUpDownTimes", function(x, min_signal, min_gap, ...) {
    standardGeneric("plotUpDownTimes")
  })
}

# One method that gets upDOwnTimes from a Stream Object
plotUpDownTimes.Stream <- function(x, min_signal=30, min_gap=60, ...) {

  upDownTimes <- getUpDownTimes(x, min_signal=min_signal, min_gap=min_gap)
  
  # NOTE:  Extra conversion step needed to guarantee display as GMT times
  GMTTimes <- as.POSIXct(upDownTimes, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  onOffs <- seq(length(GMTTimes)) %% 2 # 1, 0, 1, 0, ...
  
  # Plot the transitions with appropriate axes
  # Extend the times and onOff transitions to the full time period requested
  allTimes <- c(x@requestedStarttime, GMTTimes, x@requestedEndtime)
  allOnOffs <- c(0,onOffs,0)
  plot(allTimes, allOnOffs, type="s", xlab="GMT", ylab="", yaxt="n", ...)
  graphics::axis(2, at=c(0,1), labels=c("Off","On"), las=1, tick=TRUE)
  
  # Cover over any transition to down at the end
  # NOTE:  Seems not to work if I use GMTTimes so revert back to upDownTimes
  graphics::abline(v=x@requestedEndtime, lwd=3, col='white')
  graphics::box()
  
  # Add the title
  id <- stringr::str_sub(x@traces[[1]]@id, 1, stringr::str_length(x@traces[[1]]@id)-2)
  main <- paste("On/Off transitions for ",id)
  sensorText <- paste("(", x@traces[[1]]@Sensor, ")") 
  graphics::title(main)
  graphics::mtext(sensorText, line=0.2)
   
}

# Another method that accepts a pre-generated vector of plotUpDownTimes
plotUpDownTimes.POSIXct <- function(x, min_signal=30, min_gap=60, ...) {
  
  # NOTE:  min_signal and min_gap is not used in this function.
  # NOTE:  They are included for consistency.
  upDownTimes <- x
  
  # NOTE:  Extra step needed to guarantee display as GMT times
  GMTTimes <- as.POSIXct(upDownTimes, "%Y-%m-%dT%H:%M:%OS", tz="GMT")
  onOff <- seq(length(GMTTimes)) %% 2 # 1,0,1 0,...
  
  # TODO:  Check for and use additional plotting parameters like "main" or "sub"
  # TODO:  Use optional "starttime", "endtime" arguments to set xlim
  
  # Plot the transitions with appropriate axes
  plot(GMTTimes, onOff, type="s", xlab="GMT", ylab="", yaxt="n", ...)
  graphics::axis(2, at=c(0,1), labels=c("Off","On"), las=1, tick=TRUE)
  
  # Remove any transition to down at the end
  # NOTE:  Seems not to work if I use GMTTimes so rever back to upDownTimes
  graphics::abline(v=upDownTimes[length(upDownTimes)], lwd=2, col='white')
  graphics::box()
  
}

setMethod("plotUpDownTimes", signature(x="Stream", min_signal="numeric", min_gap="numeric"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.Stream(x, min_signal, min_gap=min_gap, ...))
setMethod("plotUpDownTimes", signature(x="Stream", min_signal="numeric", min_gap="missing"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.Stream(x, min_signal, min_gap=60, ...))
setMethod("plotUpDownTimes", signature(x="Stream", min_signal="missing", min_gap="numeric"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.Stream(x, min_signal=30, min_gap, ...))
setMethod("plotUpDownTimes", signature(x="Stream", min_signal="missing", min_gap="missing"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.Stream(x, min_signal=30, min_gap=60, ...))
setMethod("plotUpDownTimes", signature(x="POSIXct", min_signal="numeric", min_gap="numeric"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.POSIXct(x, min_signal, min_gap, ...))
setMethod("plotUpDownTimes", signature(x="POSIXct", min_signal="numeric", min_gap="missing"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.POSIXct(x, min_signal=30, min_gap, ...))
setMethod("plotUpDownTimes", signature(x="POSIXct", min_signal="missing", min_gap="numeric"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.POSIXct(x, min_signal, min_gap=60, ...))
setMethod("plotUpDownTimes", signature(x="POSIXct", min_signal="missing", min_gap="missing"), 
          function(x, min_signal, min_gap, ...) plotUpDownTimes.POSIXct(x, min_signal=30, min_gap=60, ...))

################################################################################
# Merging upDownTimes to return a new vector of datetimes that mark when
# both (either) channel is on/off
################################################################################

if (!isGeneric("mergeUpDownTimes")) {
  setGeneric("mergeUpDownTimes", function(udt1,udt2,bothOn) {
    standardGeneric("mergeUpDownTimes")
  })
}
mergeUpDownTimes.POSIXct <- function(udt1,udt2,bothOn) {
  
  # Deal with udt1 or udt2 = c()
  if (is.null(udt1)) {
    return(udt2)
  }
  if (is.null(udt2)) {
    return(udt1)
  }
  
  # Combine the times and sort them
  unsorted_times <- c(udt1,udt2)
  sort_indices <- order(unsorted_times)
  times <- unsorted_times[sort_indices]
  
  # Create a numeric vector of on/off flags with on=1, off=-1
  onOff1 <- (seq(length(udt1)) %% 2 - 0.5) * 2
  onOff2 <- (seq(length(udt2)) %% 2 - 0.5) * 2
  # Combine the flags and sort them temporally
  unsorted_onOff <- c(onOff1,onOff2)
  onOff <- unsorted_onOff[sort_indices]
  
  # NOTE:  The times and onOffs are now sorted and we will forgout about 'sort_indices'.
  # NOTE:  The indices in the next section have nothing to do with the original sort_indices
  
  # The cumulative sum now has the following meaning: 2=both on, 0=1 on, -1=both off
  cumOnOff <- cumsum(onOff)
  
  if (bothOn) {
    
    # Create up_down_times associated with "both channels on/either channel off"
    both_up_indices <- which(cumOnOff == 2)
    any_down_indices <- both_up_indices + 1
    both_upDownTimes <- sort(times[c(both_up_indices,any_down_indices)])
    return(both_upDownTimes)
    
  } else {
    
    # Create up_down_times associated with "either channel on/both channels off"
    both_down_indices <- which(cumOnOff == 0)
    # up indices start with the first time and use every time just after both down (except for the last one)
    either_up_indices <- c(1, both_down_indices[-length(both_down_indices)] + 1)
    either_upDownTimes <- sort(times[c(either_up_indices,both_down_indices)])
    return(either_upDownTimes)
    
  }
  
}
setMethod("mergeUpDownTimes", signature(udt1="POSIXct", udt2="POSIXct", bothOn="logical"), 
          function(udt1,udt2,bothOn) mergeUpDownTimes.POSIXct(udt1,udt2,bothOn))
setMethod("mergeUpDownTimes", signature(udt1="POSIXct", udt2="POSIXct", bothOn="missing"), 
          function(udt1,udt2,bothOn) mergeUpDownTimes.POSIXct(udt1,udt2,bothOn=FALSE))
setMethod("mergeUpDownTimes", signature(udt1="NULL", udt2="POSIXct", bothOn="logical"), 
          function(udt1,udt2,bothOn) mergeUpDownTimes.POSIXct(udt1,udt2,bothOn))
setMethod("mergeUpDownTimes", signature(udt1="NULL", udt2="POSIXct", bothOn="missing"), 
          function(udt1,udt2,bothOn) mergeUpDownTimes.POSIXct(udt1,udt2,bothOn=FALSE))
setMethod("mergeUpDownTimes", signature(udt1="POSIXct", udt2="NULL", bothOn="logical"), 
          function(udt1,udt2,bothOn) mergeUpDownTimes.POSIXct(udt1,udt2,bothOn))
setMethod("mergeUpDownTimes", signature(udt1="POSIXct", udt2="NULL", bothOn="missing"), 
          function(udt1,udt2,bothOn) mergeUpDownTimes.POSIXct(udt1,udt2,bothOn=FALSE))


################################################################################
# Basic plotting
################################################################################

# NOTE:  I don't seem to be able to use optional arguments of type "ANY" in the 
# NOTE:  function declaration.  Probably has to do with plot.~() being older
# NOTE:  style functions rather than S4 methods.

plot.Stream <- function(x, ...) {

  tr <- mergeTraces(x)@traces[[1]]
  
  plot(tr, ...)
  if (length(x@traces) == 1) {
    graphics::mtext(paste(length(x@traces),"trace"), side=3, line=0.5, adj=0.05)          
  } else {
    graphics::mtext(paste(length(x@traces),"traces"), side=3, line=0.5, adj=0.05)          
  }
  
}

# NOTE:  I don't seem to be able to use optional arguments of type "ANY" in the 
# NOTE:  function declaration.  Probably has to do with plot.~() being older
# NOTE:  style functions rather than S4 methods.

setMethod("plot", signature(x="Stream"), function(x, ...) plot.Stream(x, ...))
