##
##    Metric to calculate STA/LTA metric for a Stream of seismic data
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
# This metric returns the maximum value of the STA/LTA algorithm applied to this
# stream's data with the following parameters passed to the STALTA method:
#
# staSecs = 3
# ltaSecs = 30
# algorithm = "classic_RR" | "classic_LR" (default) | "classic_LR2" | "EarleAndShearer_envelope | Wong_MER"

STALTAMetric <- function(st, staSecs=3, ltaSecs=30, increment=1, algorithm="classic_LR") {
  
  starttime <- st@requestedStarttime
  endtime <- st@requestedEndtime
  
  # Make sure we're working with a single snclq
  unique_ids <- uniqueIds(st)
  if (length(unique_ids) > 1) {
    stop(paste("meanMetric: Stream has",unique_ids,"unique identifiers"))
  }
  snclq <- unique_ids[1]
  
  # Preprocessing parameters passed to STALTA
  demean <- TRUE
  detrend <- TRUE
  taper <- 0.0
  
  maxSTALTA <- 0.0
  eventTime <- starttime
  
  # Loop through all traces
  for (trace in st@traces) {
    
    # Make sure trace has enough data
    nlta <- ltaSecs * trace@stats@sampling_rate
    if (length(trace) < nlta) {
      next
    }
    
    # Get a vector of STALTA values
    stalta <- STALTA(trace, staSecs, ltaSecs, algorithm, demean, detrend, taper, increment)
    traceMaxSTALTA <- max(stalta, na.rm=TRUE)
    
    # Calculate the time at which this STALTA maximum occurred
    eventIndex <- which(stalta == traceMaxSTALTA)[1]
    traceEventTime <- trace@stats@starttime + eventIndex / trace@stats@sampling_rate    

    if (traceMaxSTALTA > maxSTALTA) {
      maxSTALTA <- traceMaxSTALTA
      eventTime <- traceEventTime
    }
    
  }
  
  # Create and return a list of Metric objects
  m1 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime,
            metricName="max_stalta", value=maxSTALTA,
            attributeName="time", attributeValueString=format(eventTime,format="%Y-%m-%dT%H:%M:%S"))
  
  return(c(m1))
  
}
