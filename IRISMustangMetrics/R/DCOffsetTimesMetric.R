##
##    Metric to determine times of DC offsets
##
##    Copyright (C) 2013  Mazama Science, Inc.
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
#

DCOffsetTimesMetric <- function(st, windowSecs=1800, incrementSecs=windowSecs/2, threshold=0.9) {
  
  # Merge traces
  tr <- mergeTraces(st,"fillNA")@traces[[1]]
  
  # Get information from the trace
  snclq <- tr@id
  starttime <- tr@stats@starttime
  endtime <- tr@stats@endtime
  
  # Bail out if we have a DC signal
  if (isDC(tr)) {
    stop(paste("SNRMetric: Trace data is a DC signal."))
  }
  
  # Calculate chunk information
  windowSamples <- windowSecs * tr@stats@sampling_rate
  incrementSamples <- incrementSecs * tr@stats@sampling_rate 
  outLength <- length(tr) / (incrementSamples)
  
  # Initialize vectors
  means <- vector("numeric",outLength) * NA
  sds <- vector("numeric",outLength) * NA
  indices <- vector("integer",outLength) * NA
  
  # Loop through all chunks
  for (i in seq(outLength)) {
    lo <- (i-1) * incrementSamples + 1
    hi <- lo + windowSamples
    means[i] <- mean(tr@data[lo:hi],na.rm=TRUE)
    sds[i] <- sd(tr@data[lo:hi],na.rm=TRUE)
    indices[i] <- as.integer(lo)
  }  

  # NOTE:  diff() returns a vector that is one less than the incoming vector.
  # NOTE:  Thus, diff(mean)[1] = mean[2]-mean[1].
  # NOTE:  indices[N] referes to the location at the beginning of mean[N] so we
  # NOTE:  really want to associate indices[N] with diff(mean[N+1])
  
  metric <- abs(diff(means)) / max(mean(sds,na.rm=TRUE),1)
  jumps <- which(metric > threshold) + 1 # because metric is one shorter than indices
  dcOffsetTimes <- starttime + indices[jumps]/tr@stats@sampling_rate
  
  # Create and return a MultipleTimeValue metric
  m1 <- new("MultipleTimeValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="dc_offset_times", values=dcOffsetTimes)

  return(c(m1))
  
}
