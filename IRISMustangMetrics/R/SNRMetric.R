##
##    Metric to calculate a signal/noise ratio comparing the rmsVariance (power)
##    by dividing the incoming stream into two equal parts.  It is assumed that
##    the stream starttime and endtime are equally spaced about the onset of
##    of a seismic event.
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
#

SNRMetric <- function(st, algorithm="splitWindow", windowSecs=600) {
    
  # Extract the trace and associated information
  tr <- st@traces[[1]]
  snclq <- tr@id
  starttime <- tr@stats@starttime
  endtime <- tr@stats@endtime
  
  # Chck for gaps
  if (length(st@traces) > 1) {    
    stop(paste("skipping", snclq, "because it has gaps"))
  }
     
  # Bail out if we have don't have enough data
  if (as.numeric(difftime(endtime, starttime, units="sec")) < windowSecs) {
    stop(paste("SNRMetric: Data do not fill the window."))
  }

  # Bail out if we have a DC signal
  if (isDC(tr)) {
    stop(paste("SNRMetric: Trace data is a DC signal."))
  }
  
  if (algorithm == "splitWindow") {
    
    # NOTE:  In this case, the window is assumed to be centered about
    # NOTE:  an event, perhaps determined withthe IRIS DMC "event" webservice
    # NOTE:  The triggerOnset is just the midpoint of the trace.

    to <- starttime + as.numeric(difftime(endtime, starttime, units="sec")) / 2
    
    trN <- slice(tr, to-windowSecs/2, to)
    trS <- slice(tr, to, to+windowSecs/2)
    snr <- rmsVariance(trS) / rmsVariance(trN)
    
  } else if (algorithm == "staltaTrigger") {
    
    # Demean and detrend the data first
    tr <- DDT(tr, TRUE, TRUE, 0)

    # Find the P-wave onset with "classic_LR"
    staSecs <- 3
    ltaSecs <- 30 
    # Make sure that there are at least two points in the STA window
    if (staSecs * tr@stats@sampling_rate < 1) {
      staSecs <- 2 / tr@stats@sampling_rate
      ltaSecs <- staSecs * 5
    }
    picker <- STALTA(tr,staSecs,ltaSecs,"classic_LR")
    to <- triggerOnset(tr,picker)
    
    trN <- slice(tr, to-windowSecs/2, to)
    trS <- slice(tr, to, to+windowSecs/2)
    snr <- rmsVariance(trS) / rmsVariance(trN)
    
  }
  
  # Create and return a list of Metric objects
  m1 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="sample_snr", value=snr)
  
  return(c(m1))
  
}
