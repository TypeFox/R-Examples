##
##    S4 classes for handling seismic traces
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
# R classes for a Trace obj.
#
# This is a port of some of the functionality found in obspy.core.trace.py 
# 
#   http://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.html
################################################################################

################################################################################
# Class TraceHeader
#
# A container for additional header informaion for each Trace object.
#
# A TraceHeader object contains all header information (aka metadata)
# associated with a Trace object.
#
# Documentation for attributes copied verbatim from obspy.core.trace.py:
# 
#    ``sampling_rate`` : float, optional
#        Sampling rate in hertz (default value is 1.0).
#    ``delta`` : float, optional
#        Sample distance in seconds (default value is 1.0).
#    ``calib`` : float, optional
#        Calibration factor (default value is 1.0).
#    ``npts`` : int, optional
#        Number of sample points (default value is 0, which implies that no data
#        is present).
#    ``network`` : string, optional
#        Network code (default is an empty string).
#    ``location`` : string, optional
#        Location code (default is an empty string).
#    ``station`` : string, optional
#        Station code (default is an empty string).
#    ``channel`` : string, optional
#        Channel code (default is an empty string).
#    ``starttime`` : :class:`~obspy.core.utcdatetime.UTCDateTime`, optional
#        Date and time of the first data sample given in UTC (default value is
#        "1970-01-01T00:00:00.0Z").
#    ``endtime`` : :class:`~obspy.core.utcdatetime.UTCDateTime`, optional
#        Date and time of the last data sample given in UTC
#        (default value is "1970-01-01T00:00:00.0Z").
#
################################################################################

setClass("TraceHeader", 
  # typed slots (aka attributes) for class TraceHeader
  representation(sampling_rate = "numeric",
                 delta = "numeric",
                 calib = "numeric",
                 npts = "integer",
                 network = "character",
                 location = "character",
                 station = "character",
                 channel = "character",
                 quality = "character",
                 starttime = "POSIXct",
                 endtime = "POSIXct",
                 processing = "character"),
  # default values for slots
  prototype(sampling_rate = 1.0,
            delta = 1.0,
            calib = 1.0,
            npts = as.integer(0),
            network = "",
            location = "",
            station = "",
            channel = "",
            quality = "X",
            starttime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"),
            endtime = as.POSIXct("1900-01-01T00:00:00",format="%Y-%m-%dT%H:%M:%OS", tz="GMT"))
)

################################################################################
# TODO:  Improve notes on TraceHeader initialization
#
# A TraceHeader can be initialized with a headerLine as returned by the IRIS DMC timeseries webservice:
#
# TIMESERIES LD_POTS__HHZ_M, 351 samples, 100.503 sps, 2012-01-29T00:00:00.006000, SLIST, INTEGER, COUNTS
#
# or with named elements from a headerList.  Items in the headerList will override items determined
# by parsing the headerLine.
################################################################################

# TODO:  I am currently using "Trace" to mean what libmseed calls a "Trace Segment".
# TODO:  Each libmseed "Trace" has both "numsamples" and "samplecnt".
# TODO:  Need to figure out what the difference is and when to bail with zero samples in the SNCL.
# TODO:  Zero samples in a "Trace Segment" is not necessariliy a deal breaker.
setMethod("initialize", "TraceHeader",
  function(.Object,
           headerList=list(),
           headerLine=NULL,
           ...) {
    
    # Begin by parsing any headerLine that was passed in.
    if (!is.null(headerLine)) {
      headerValues <- strsplit(headerLine,", ")[[1]]
      sncl <- strsplit(headerValues[1]," ")[[1]][2]
      snclValues <- strsplit(sncl,"_")[[1]]
      .Object@network <- snclValues[1]
      .Object@station <- snclValues[2]
      .Object@location <- snclValues[3]
      .Object@channel <- snclValues[4]
      .Object@npts <- as.integer(strsplit(headerValues[2]," ")[[1]][1])
      .Object@sampling_rate <- as.numeric(strsplit(headerValues[3]," ")[[1]][1])
      .Object@starttime <- as.POSIXct(headerValues[4], format="%Y-%m-%dT%H:%M:%OS", tz="GMT") # header uses 'T' format
    }
    
    # override with anything found in the headerList
    for (key in names(headerList)) {
      slot(.Object,key) = headerList[[key]]
    }
  
    # guarantee npts and sampling_rate are non-zero numbers
    if (is.na(.Object@npts)) {
      err_msg <- paste("initialize.TraceHeader: npts=", .Object@npts, sep="")
      if (!is.null(headerLine)) {
        err_msg <- paste(err_msg, "\nheaderLine =\"", headerLine, "\"", sep="")
      }
      stop(err_msg)
    }
    if (is.na(.Object@sampling_rate) || (.Object@sampling_rate == 0)) {
      err_msg <- paste("initialize.TraceHeader: sampling_rate=", .Object@sampling_rate, sep="")
      if (!is.null(headerLine)) {
        err_msg <- paste(err_msg, "\nheaderLine =\"", headerLine, "\"", sep="")
      }
      stop(err_msg)
    }    
    
    # calculate derived values    
    delta = 1.0 / .Object@sampling_rate
    if (!is.finite(delta)) {
      delta = 0
    }
    .Object@delta = delta

    if (.Object@npts == 0) {
      delta = 0
    } else {
      delta = (.Object@npts-1) / .Object@sampling_rate
      if (!is.finite(delta))
        delta = 0
    }
    .Object@endtime = .Object@starttime + delta

    # return the newly created .Object
    return(.Object)
  }
)


# show method is called by show() and print() ----------------------------------

show.TraceHeader <- function(object) {
  cat ("Seismic Trace TraceHeader \n")
  cat (" Network:       " , object@network, "\n")
  cat (" Station:       " , object@station, "\n")
  cat (" Location:      " , object@location, "\n")
  cat (" Channel:       " , object@channel, "\n")
  cat (" Quality:       " , object@quality, "\n")
  cat (" calib:         " , object@calib, "\n")
  cat (" npts:          " , object@npts, "\n")
  cat (" sampling rate: " , object@sampling_rate, "\n")
  cat (" delta:         " , object@delta, "\n")
  cat (" starttime:     " , format(object@starttime), "\n")
  cat (" endtime:       " , format(object@endtime), "\n")
  cat (" processing:    " , paste(object@processing,collapse="; "), "\n")
}
# NOTE:  method signature must match generic signature for 'show' with argument: 'object'
setMethod("show", signature(object="TraceHeader"), function(object) show.TraceHeader(object))


# as.headerLine prints TraceHeader out as a single line ------------------------ 
#
# TIMESERIES LD_POTS__HHZ_M, 351 samples, 100.503 sps, 2012-01-29T00:00:00.006000, SLIST, INTEGER, COUNTS

if (!isGeneric("as.headerLine")) {
  setGeneric("as.headerLine", function(obj) {
    standardGeneric("as.headerLine")
  })
}
as.headerLine.TraceHeader <- function(obj) {
  cat("TIMESERIES ",obj@network,"_",obj@station,"_",obj@location,"_",obj@channel,", ", sep="")
  cat(obj@npts," samples, ",obj@sampling_rate," sps, ",obj@starttime,", SLIST, INTEGER, COUNTS", sep="")
}
setMethod("as.headerLine", signature(obj="TraceHeader"), function(obj) as.headerLine.TraceHeader(obj))

# TODO:  Setter and getter methods could be written as described in:
# TODO:    www.stat.auckland.ac.nz/S-Workshop/Gentleman/S4Objects.pdf


################################################################################
# Class Trace
#
# An object containing data of a continuous series, such as a seismic trace.
#
# Documentation for attributes copied verbatim from obspy.core.trace.py
# 
#    :type data: :class:`~numpy.ndarray` or :class:`~numpy.ma.MaskedArray`
#    :param data: Array of data samples
#    :type header: dict or :class:`~obspy.core.trace.Stats`
#    :param header: Dictionary containing header fields
#
#    :var id: A SEED compatible identifier of the trace.
#    :var stats: A container :class:`~obspy.core.trace.Stats` for additional
#        header information of the trace.
#    :var data: Data samples in a :class:`~numpy.ndarray` or
#        :class:`~numpy.ma.MaskedArray`
#
################################################################################

setClass("Trace", 
  # typed slots (aka attributes) for class Trace
  representation(id = "character",
                 stats = "TraceHeader",
                 Sensor = "character",
                 InstrumentSensitivity = "numeric",
                 InputUnits = "character",
                 data = "numeric"),
  # default values for slots
  prototype(id = "",
            stats = new("TraceHeader"),
            Sensor = "",
            InstrumentSensitivity = 1.0,
            InputUnits = "",
            data = c(0))
)

# initialze method
setMethod("initialize", "Trace",
  function(.Object,
           id="",
           stats=new("TraceHeader"),
           Sensor = "",
           InstrumentSensitivity = 1.0,
           InputUnits = "",
           data=c(0),
           ...) {
    
    .Object@id <- id
    .Object@stats <- stats
    .Object@Sensor <- Sensor
    .Object@InstrumentSensitivity <- InstrumentSensitivity
    .Object@InputUnits <- InputUnits
    .Object@data <- data

    if (.Object@id == "") {
      stats <- .Object@stats
      .Object@id <- paste(stats@network,stats@station,stats@location,stats@channel,stats@quality, sep=".")
    }

    return(.Object)
  }
)


################################################################################
# Simple methods working on the data from a Trace object
################################################################################

# NOTE:  The method signatures must match those already defined by R core functions
# NOTE:  and R standard generics.  Use "?function" to see the documentation or just
# NOTE:  "function" to see the implementation and signature.
# NOTE:
# NOTE:  When creating a new method to extend a method that has "...", always
# NOTE:  provide the S4 object as the first argument, then the "..." and then 
# NOTE:  any final arguments.


# as.vector --------------------------------------------------------------------

as.vector.Trace <- function(x, mode="any") {
  return( x@data )
}
setMethod("as.vector", signature(x="Trace"), function(x) as.vector.Trace(x))


# isDC -------------------------------------------------------------------------

if (!isGeneric("isDC")) {
  setGeneric("isDC", function(x) {
    standardGeneric("isDC")
  })
} 
isDC.Trace <- function(x) {
  x <- DDT(x, demean=TRUE, detrend=FALSE, taper=0)
  return( all(x@data==0, na.rm=TRUE) )
}
setMethod("isDC", signature(x="Trace"), function(x) isDC.Trace(x))


# Length -----------------------------------------------------------------------

length.Trace <- function(x) {
  return( length(x@data) )
}
setMethod("length", signature(x="Trace"), function(x) length.Trace(x))


# Maximum ----------------------------------------------------------------------

max.Trace <- function(x, ..., na.rm=TRUE) {
  return( max(x@data, na.rm=na.rm) )
}
setMethod("max", signature(x="Trace"), function(x, ...) max.Trace(x, ...))


# Mean -------------------------------------------------------------------------

mean.Trace <- function(x, ...) {
  return( mean(x@data, ...) )
}
setMethod("mean", signature(x="Trace"), function(x, ...) mean.Trace(x, ...))
# TODO:  How to get mean to use na.rm=TRUE as the default?


# Median -----------------------------------------------------------------------

median.Trace <- function(x, na.rm) {
  return( median(x@data, na.rm=na.rm) )
}
# NOTE:  method signature must match generic signature for function 'median' with arguments: 'x', 'na.rm'
setMethod("median", signature(x="Trace", na.rm="logical"), function(x, na.rm) median.Trace(x, na.rm=na.rm))
setMethod("median", signature(x="Trace", na.rm="missing"), function(x, na.rm) median.Trace(x, na.rm=TRUE))


# Minimum ----------------------------------------------------------------------

min.Trace <- function(x, ..., na.rm=TRUE) {
  return( min(x@data, na.rm=na.rm) )
}
setMethod("min", signature(x="Trace"), function(x, ...) min.Trace(x, ...))


# Standard Deviation -----------------------------------------------------------

sd.Trace <- function(x, na.rm) {
  return( sd(x@data, na.rm=na.rm) )
}
# NOTE:  method signature must match generic signature for for function 'sd' with arguments: 'x', 'na.rm'
setMethod("sd", signature(x="Trace", na.rm="logical"), function(x, na.rm) sd.Trace(x, na.rm=na.rm))
setMethod("sd", signature(x="Trace", na.rm="missing"), function(x, na.rm) sd.Trace(x, na.rm=TRUE))


# Root Mean Square ----------------------------------------------------

if (!isGeneric("rms")) {
  setGeneric("rms", function(x, na.rm) {
    standardGeneric("rms")
  })
} 
rms.Trace <- function(x, na.rm) {
  return( sqrt( mean((x@data)^2, na.rm=na.rm) ) )
}
setMethod("rms", signature(x="Trace", na.rm="logical"), function(x, na.rm) rms.Trace(x, na.rm=na.rm))
setMethod("rms", signature(x="Trace", na.rm="missing"), function(x, na.rm) rms.Trace(x, na.rm=TRUE))


# Root Mean Square Variance ----------------------------------------------------

if (!isGeneric("rmsVariance")) {
  setGeneric("rmsVariance", function(x, na.rm) {
    standardGeneric("rmsVariance")
  })
} 
rmsVariance.Trace <- function(x, na.rm) {
  mean <- mean(x, na.rm=na.rm)
  n <- length(x)
  return( sqrt( sum( (x@data-mean)^2, na.rm=na.rm ) / (n) ) )
}
setMethod("rmsVariance", signature(x="Trace", na.rm="logical"), function(x, na.rm) rmsVariance.Trace(x, na.rm=na.rm))
setMethod("rmsVariance", signature(x="Trace", na.rm="missing"), function(x, na.rm) rmsVariance.Trace(x, na.rm=TRUE))


################################################################################
# Methods that return a new, modified Trace
################################################################################

# Multiply by a constant -------------------------------------------------

if (!isGeneric("multiplyBy")) {
  setGeneric("multiplyBy", function(x, y) {
    standardGeneric("multiplyBy")
  })
} 

multiplyBy.Trace <- function(x, y) {
  id <- x@id
  stats <- x@stats
  Sensor <- x@Sensor
  InstrumentSensitivity <- x@InstrumentSensitivity
  InputUnits <- x@InputUnits
  data <- x@data * y
  stats@processing <- append(stats@processing,paste("multiply by",y))
  
  return( new("Trace", id, stats, Sensor, InstrumentSensitivity, InputUnits, data) )
}

setMethod("multiplyBy", signature("Trace", y="numeric"), function(x, y) multiplyBy.Trace(x, y=y))


# Apply typical seismic DDT ----------------------------

if (!isGeneric("DDT")) {
  setGeneric("DDT", function(x, demean, detrend, taper) {
    standardGeneric("DDT")
  })
} 

DDT.Trace <- function(x, demean, detrend, taper) {
  
  # TODO:  When traces have processingInfo we can check the processingInfo
  # TODO:  and skip any steps that have already been done.
  
  id <- x@id
  stats <- x@stats
  Sensor <- x@Sensor
  InstrumentSensitivity <- x@InstrumentSensitivity
  InputUnits <- x@InputUnits
  
  # NOTE:  Use the pracma::detrend() function.
  # NOTE:  Linear detrending also removes the mean.
  
  # demean only
  if (demean && !detrend) {
    data <- x@data - mean(x@data, na.rm=TRUE)
    stats@processing <- append(stats@processing,"demean")
  } else {
    data <- x@data
  }
  
  # detrend
  if (detrend) {
    data <- as.numeric(pracma::detrend(x@data, tt='linear'))
    stats@processing <- append(stats@processing,"detrend")
  }
  
  # cosine taper
  if (taper > 0) {
    data <- stats::spec.taper(data, p=taper)
    stats@processing <- append(stats@processing,"cosine taper")
  }
  

  return( new("Trace", id, stats, Sensor, InstrumentSensitivity, InputUnits, data) )
}

# All parameters specified
setMethod("DDT", signature("Trace", demean="logical", detrend="logical", taper="numeric"),
          function(x, demean, detrend, taper) DDT.Trace(x, demean, detrend, taper))
# No parameters
setMethod("DDT", signature("Trace", demean="missing", detrend="missing", taper="missing"),
          function(x, demean, detrend, taper) DDT.Trace(x, TRUE, TRUE, 0.1))


# Apply Butterworth filters ------------

if (!isGeneric("butterworth")) {
  setGeneric("butterworth", function(x, n, low, high, type) {
    standardGeneric("butterworth")
  })
} 

butterworth.Trace <- function(x, n, low, high, type) {
  
  # TODO:  When traces have processingInfo we can check the processingInfo
  # TODO:  and skip any steps that have already been done.
  
  data <- x@data
  id <- x@id
  stats <- x@stats
  Sensor <- x@Sensor
  InstrumentSensitivity <- x@InstrumentSensitivity
  InputUnits <- x@InputUnits
  
  # If neither demean nor detrend has ever been applied, demean the data
  if ( !any(stringr::str_detect(x@stats@processing,"demean")) &&
       !any(stringr::str_detect(x@stats@processing,"detrend")) ) {
    data <- data - mean(data, na.rm=TRUE)
    stats@processing <- append(stats@processing,"demean")
  }
  
  
#   filt <- signal::butter(2,c(0.02/(trZ@stats@sampling_rate/2),0.04/(trZ@stats@sampling_rate/2)),type="pass")
#   tsZ <- stats::ts(trZ@data,frequency=trZ@stats@sampling_rate)
#   ts2 <- stats::ts(tr2@data,frequency=tr2@stats@sampling_rate)
#   ts1 <- stats::ts(tr1@data,frequency=tr1@stats@sampling_rate)
#   
#   trZ_f <- trZ; trZ_f@data <- as.numeric(signal::filter(filt,tsZ))
#   
  
  # Create normalization factor
  norm <- stats@sampling_rate / 2
  
  # Create Butterworth filter
  if (type == 'low') {
    bf <- signal::butter(n, low/norm, type='low')
  } else if (type == 'high') {
    bf <- signal::butter(n, high/norm, type='high')
  } else {
    bf <- signal::butter(n, c(low/norm,high/norm), type=type)
  }
  
  # Apply filter
  ts <- stats::ts(data, frequency=stats@sampling_rate)
  ts_filtered <- signal::filter(bf, ts)
  data <- as.numeric(ts_filtered)
  stats@processing <- append(stats@processing,paste0("Butterworth filtered with (",n,",",low,",",high,",",type,")"))
  
  return( new("Trace", id, stats, Sensor, InstrumentSensitivity, InputUnits, data) )
}

# All parameters specified
setMethod("butterworth", signature("Trace", n="numeric", low="numeric", high="numeric", type="character"),
          function(x, n, low, high, type) butterworth.Trace(x, n, low, high, type))
# type missing -- band pass
setMethod("butterworth", signature("Trace", n="numeric", low="numeric", high="numeric", type="missing"),
          function(x, n, low, high, type) butterworth.Trace(x, n, low, high, "pass"))
# hi and type missing -- low pass
setMethod("butterworth", signature("Trace", n="numeric", low="numeric", high="missing", type="missing"),
          function(x, n, low, high, type) butterworth.Trace(x, n, low, NULL, "high"))
# lo and type missing -- high pass
setMethod("butterworth", signature("Trace", n="numeric", low="missing", high="numeric", type="missing"),
          function(x, n, low, high, type) butterworth.Trace(x, n, NULL, high, "low"))


# Slice out a portion of a Trace ---------------------------

# NOTE:  Unlike the ObsPy method, this method does no padding.  The returned 
# NOTE:  Trace will always be a subset of the original Trace.  Slicing is rounded
# NOTE:  to the nearest second.

if (!isGeneric("slice")) {
  setGeneric("slice", function(x, starttime, endtime) {
    standardGeneric("slice")
  })
} 
slice.Trace <- function(x, starttime, endtime) {
  
  # id and stats will be part of the returned Trace
  id <- x@id
  stats <- x@stats
  Sensor <- x@Sensor
  InstrumentSensitivity <- x@InstrumentSensitivity
  InputUnits <- x@InputUnits
  data <- x@data

  sliced <- FALSE
  
  start_index <- 1
  end_index <- length(data)
  
  # Sanity check
  if (starttime >= endtime) {
    stop(paste("slice.Trace: requested starttime \"", starttime, "\" >= requested endtime \"", endtime, "\""))    
  }
  if (starttime >= stats@endtime) {
    stop(paste("slice.Trace: requested starttime \"", starttime, "\" >= Trace endtime \"", stats@endtime, "\""))
  }
  if (endtime <= stats@starttime) {
    stop(paste("slice.Trace: requested endtime \"", endtime, "\" <= Trace starttime \"", stats@starttime, "\""))
  }
  
  # Modify the startime and startIndex if needed
  if (starttime > stats@starttime ) {    
    # NOTE:  Use difftime() instead of just subtracting to guarantee that units are "secs"
    offset_secs <- as.numeric(round( difftime(starttime, stats@starttime, units="secs") ))
    start_index <- start_index + as.integer(floor(offset_secs * stats@sampling_rate))
    stats@starttime <- starttime
    sliced <- TRUE
  } # else leave the start as is
  
  # Modify the endtime and end_index if needed
  if (endtime < stats@endtime) {    
    # NOTE:  Use difftime() instead of just subtracting to guarantee that units are "secs"
    offset_secs <- as.numeric(round( difftime(stats@endtime, endtime, units="secs") ))
    end_index <- end_index - as.integer(floor(offset_secs * stats@sampling_rate))
    stats@endtime <- endtime       
    sliced <- TRUE
  } # else leave the end as is
  
  if (sliced) {
    stats@npts <- as.integer(end_index - start_index + 1)
    stats@processing <- append(stats@processing,"slice")
    return( new("Trace", id, stats, Sensor, InstrumentSensitivity, InputUnits, data=data[start_index:end_index]) )    
  } else {
    return(x)
  }
  
}
setMethod("slice", signature(x="Trace", starttime="POSIXct", endtime="POSIXct"),
          function(x, starttime, endtime) slice.Trace(x, starttime=starttime, endtime=endtime))



################################################################################
# Spectral Methods
################################################################################

# Hilbert transform -------------------------

if (!isGeneric("hilbert")) {
  setGeneric("hilbert", function(x) {
    standardGeneric("hilbert")
  })
} 

# Code copied from seewave package hilbert() method
#   http://cran.r-project.org/web/packages/seewave/index.html
hilbert.Trace <- function(x) {
  
  # Always detrend, demean, cosine taper the data
  x <- DDT(x)
  hfft <- hilbertFFT(x@data)
  # The Hilbert tansform is the Imaginary part of the transformed data
  x@data <- Im(hfft)
  x@stats@processing <- append(x@stats@processing,"Hilbert transform")
  
  return(x)
}

setMethod("hilbert", signature(x="Trace"),
          function(x) hilbert.Trace(x))



# Envelope using Hilbert transform -------------------------

if (!isGeneric("envelope")) {
  setGeneric("envelope", function(x) {
    standardGeneric("envelope")
  })
} 

# Code copied from seewave package env() and hilbert() methods
#   http://cran.r-project.org/web/packages/seewave/index.html
envelope.Trace <- function(x) {
  
  # Always detrend, demean, cosine taper the data
  x <- DDT(x)
  hfft <- hilbertFFT(x@data)
  # The envelope is the Modulus of the transformed data
  x@data <- Mod(hfft)
  x@stats@processing <- append(x@stats@processing,"envelope")
  
  return(x)
}

setMethod("envelope", signature(x="Trace"),
          function(x) envelope.Trace(x))


################################################################################
# Methods for event picking
################################################################################

# Multi-purpose SLA/TLA trigger algorithm --------------------------------------

# NOTE:  http://en.wikipedia.org/wiki/First_break_picking
# NOTE:  
# NOTE:  "A Comparison of Select Trigger Algorithms for Automated Global Seismic Phase and Event Detection"
# NOTE:  http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.116.245&rep=rep1&type=pdf
# NOTE:
# NOTE:  Wong:
# NOTE:  "Automatic time-picking of first arrivals on noisy microseismic data"
# NOTE:  http://www.crewes.org/ForOurSponsors/ConferenceAbstracts/2009/CSEG/Wong_CSEG_2009.pdf
# NOTE:
# NOTE:  GEOPHYSICS, VOL. 75, NO. 4 (JULY-AUGUST 2010); P. V67-V76
# NOTE:  Automatic first-breaks picking: New strategies and algorithms"
# NOTE:  http://www.fcaglp.unlp.edu.ar/~velis/papers/PickingGeop10.pdf
# NOTE:
# NOTE:  "Adaptive microseismic event detection and automatic time picking"
# NOTE:  http://www.cspg.org/documents/Conventions/Archives/Annual/2012/279_GC2012_Adaptive_Microseismic_Event_Detection.pdf
# NOTE:
# NOTE:  Earle and Shearer:
# NOTE:  "Characterization of Global Seismograms Using an Automatic-Picking Algorithm"
# NOTE:  Bulletin of the Seismological Society of America, Vol. 84, No. 2, pp. 366-376, April 1994

if (!isGeneric("STALTA")) {
  setGeneric("STALTA", function(x, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment) {
    standardGeneric("STALTA")
  })
} 

STALTA.Trace <- function(x, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment) {
  # Calculate constants
  n_sta <- staSecs * x@stats@sampling_rate
  n_lta <- ltaSecs * x@stats@sampling_rate
  
  if (n_sta < 1) {
    stop(paste("STALTA.Trace: STA window of",staSecs,"secs is too short for trace sampling rate of", x@stats@sampling_rate,"Hz."))
  }

  if (demean | detrend) {
    x <- DDT(x, demean, detrend, taper)   
  }
  len <- length(x)
  
  # Stop if there is insufficient data
  if ( len < n_lta ) {
    stop(paste("STALTA.Trace: insufficient data."))
  }
  
  # DC signals will have all zeroes after demeaning.
  # Return all zeroes in this case.
  if ( all(x@data == 0, na.rm=TRUE) ) {
    return(x@data)
  }
  
  # Apply the algorithm over the entire data series
    
  if (algorithm == "classic_RR") {
    # "Automatic time-picking of first arrivals on noisy microseismic data"
    # NOTE:  In this reference, the LTA and STA windows use the same right edge (index is in bth windows)
    if (increment != 1) {
      stop(paste("STALTA.Trace: algorithm EarleAndShearer_envelope requires increment=1.",sep=""))
    }
    sta <- seismicRoll::roll_mean(x@data^2, n_sta, increment, align="right")
    lta <- seismicRoll::roll_mean(x@data^2, n_lta, increment, align="right")
    picker <- sta/lta
  } else if (algorithm == "classic_LR") {
    # NOTE:  This is the same algorithm but with left-right alignment (index is in both windows)
    picker <- seismicRoll::roll_stalta(x@data^2, n_sta, n_lta, increment)
  } else if (algorithm == "EarleAndShearer_envelope") {
    if (increment != 1) {
      stop(paste("STALTA.Trace: algorithm EarleAndShearer_envelope requires increment=1.",sep=""))
    }
    data <- envelope(x)@data
    sta <- seismicRoll::roll_mean(data, n_sta, increment, align="left")
    lta <- seismicRoll::roll_mean(data, n_lta, increment, align="right")    
    picker <- sta/lta
  } else {
    stop(paste("STALTA.Trace: algorithm=\"",algorithm,"\" not recognized.",sep=""))    
  }
  
  return( picker )
}

# All parameters specified
setMethod("STALTA", signature(x="Trace", staSecs="numeric", ltaSecs="numeric",
                              algorithm="character", demean="logical", detrend="logical", taper="numeric", increment="numeric"),
          function(x, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment) 
            STALTA.Trace(x, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment))
# demean, detrend, taper and increment missing
setMethod("STALTA", signature(x="Trace", staSecs="numeric", ltaSecs="numeric",
                              algorithm="character", demean="missing", detrend="missing", taper="missing", increment="missing"),
          function(x, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment) 
            STALTA.Trace(x, staSecs, ltaSecs, algorithm, TRUE, TRUE, 0.0, 1))
# Only lstaSecs and ltaSecs
setMethod("STALTA", signature(x="Trace", staSecs="numeric", ltaSecs="numeric",
                              algorithm="missing", demean="missing", detrend="missing", taper="missing", increment="missing"),
          function(x, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment) 
            STALTA.Trace(x, staSecs, ltaSecs, "classic_LR", TRUE, TRUE, 0.0, 1))
# All parameters missing
setMethod("STALTA", signature(x="Trace", staSecs="missing", ltaSecs="missing",
                              algorithm="missing", demean="missing", detrend="missing", taper="missing", increment="missing"),
          function(x, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment) 
            STALTA.Trace(x, 3, 30, "classic_LR", TRUE, TRUE, 0.0, 1))


# Simple triggering of p-wave onset --------------------------------------------

# NOTE:  Method to return the first timepoint at which STALTA is above a threshold
#
# NOTE:  This method only returns the first timepoint or NA.  It DOES NOT implement
# NOTE:  the full functionality found in the ObsPy method of the same name.
#
# NOTE:  http://docs.obspy.org/packages/autogen/obspy.signal.trigger.triggerOnset.html

if (!isGeneric("triggerOnset")) {
  setGeneric("triggerOnset", function(x, picker, threshold, index) {
    standardGeneric("triggerOnset")
  })
} 

triggerOnset.Trace <- function(x, picker, threshold, index) {
  
  if (missing(threshold)) {
    threshold <- stats::quantile(picker,0.99999,na.rm=TRUE)
  }
  
  if (missing(index)) {
    index <- FALSE
  }
  
  # At which index is the picker above the threshold?
  eventIndex <- which(picker > threshold)[1]
  
  # Return just the index if desired
  if (index) {
    return(eventIndex)
  }
  
  # Otherwise calculate the time at which that index occurs
  if (is.na(eventIndex)) {
    eventTime <- NA
  } else {
    eventTime <- x@stats@starttime + eventIndex / x@stats@sampling_rate    
  }

  return(eventTime)
  
}

setMethod("triggerOnset", signature(x="Trace", picker="numeric", threshold="ANY", index="ANY"),
          function(x, picker, threshold, index)
            triggerOnset.Trace(x, picker, threshold, index))


# Select the event window around p-wave onset ----------------------------------

if (!isGeneric("eventWindow")) {
  setGeneric("eventWindow", function(x, picker, threshold, windowSecs) {
    standardGeneric("eventWindow")
  })
} 

eventWindow.Trace <- function(x, picker, threshold, windowSecs) {
  
  eventOnset <- triggerOnset(x, picker, threshold)
  windowStarttime <- eventOnset - windowSecs/2
  windowEndtime <- eventOnset + windowSecs/2
  
  windowTrace <- slice(x, windowStarttime, windowEndtime)
  
  return(windowTrace)
  
}

setMethod("eventWindow", signature(x="Trace", picker="numeric", threshold="numeric", windowSecs="numeric"),
          function(x, picker, threshold, windowSecs) 
            eventWindow.Trace(x, picker, threshold, windowSecs))
setMethod("eventWindow", signature(x="Trace", picker="numeric", threshold="numeric", windowSecs="missing"),
          function(x, picker, threshold, windowSecs) 
            eventWindow.Trace(x, picker, threshold, windowSecs=3600))
setMethod("eventWindow", signature(x="Trace", picker="numeric", threshold="missing", windowSecs="numeric"),
          function(x, picker, threshold, windowSecs) 
            eventWindow.Trace(x, picker, threshold=stats::quantile(picker,0.999,na.rm=TRUE), windowSecs))
setMethod("eventWindow", signature(x="Trace", picker="numeric", threshold="missing", windowSecs="missing"),
          function(x, picker, threshold, windowSecs) 
            eventWindow.Trace(x, picker, threshold=stats::quantile(picker,0.999,na.rm=TRUE), windowSecs=3600))


################################################################################
# Basic plotting
################################################################################

# NOTE:  I don't seem to be able to use optional arguments of type "ANY" in the 
# NOTE:  function declaration.  Probably has to do with plot.~() being older
# NOTE:  style functions rather than S4 methods.

plot.Trace <- function(x, starttime=x@stats@starttime, endtime=x@stats@endtime, 
                       subsampling=NULL, add=FALSE, ...) {
  
  # Create subsampled indices if necessary
  max_size <- 10000
  if (is.null(subsampling)) {
    if (length(x) <= max_size) {
      subsampling <- 1
    } else {
      subsampling <- as.integer(length(x)/(max_size/2))
    }
  }
  indices <- seq(1,length(x),subsampling)
  
  # Create labels
  
  # remove ".M" or other quality factor designation as people don't expect it
  id <- stringr::str_sub(x@id, 1, stringr::str_length(x@id)-2)
  main <- paste("Seismic Trace for ",id)
  sensorText <- paste("(", x@Sensor, ")")
  # Create array of times
  times <- seq(from=x@stats@starttime, to=x@stats@endtime, length.out=length(x@data))

  if (length(x) == length(indices)) {
    xlab <- paste("(",length(x), "points )")    
  } else {    
    xlab <- paste("( ",length(x), " points, subsampling=", subsampling, " for this plot. )", sep="")
  }
###  ylab <- paste(x@InputUnits)
  ylab <- "raw"
  
#   # Set up the time range to plot
#   if (missing(starttime)) {
#     starttime <- x@stats@starttime
#   }
#   if (missing(endtime)) {
#     endtime <- x@stats@endtime
#   }
  
  # Plot
  plot(x@data[indices] ~ times[indices], type='l', main="", #xaxt="n",
       xlim=c(starttime, endtime),
       xlab=xlab, ylab=ylab, ...)
  # x-axis
#   graphics::axis(1, at=c(1,length(indices)), labels=c("",""))
#   graphics::mtext(starttime, side=1, line=0.7, at=1, adj=0)
#   graphics::mtext(endtime, side=1, line=0.7, at=length(indices), adj=1)
  # title
  graphics::title(main)
  graphics::mtext(sensorText, line=0.2)
  
}

# NOTE:  I don't seem to be able to use optional arguments of type "ANY" in the 
# NOTE:  function declaration.  Probably has to do with plot.~() being older
# NOTE:  style functions rather than S4 methods.

setMethod("plot", signature(x="Trace"), function(x, ...) plot.Trace(x, ...))


################################################################################
# END
################################################################################



