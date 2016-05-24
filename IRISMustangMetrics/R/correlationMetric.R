##
##    Metric to calculate the correlation between two streams of seismic data
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

correlationMetric <- function(st1, st2) {
  
  # Get the first traces 
  tr1 <- st1@traces[[1]]
  tr2 <- st2@traces[[1]]
  
  starttime <- st1@requestedStarttime
  endtime <- st1@requestedEndtime

  # Sanity check stream lengths
  if (st2@requestedStarttime != starttime || st2@requestedEndtime != endtime) {
      stop(paste("correlationMetric: Incompatible starttimes or endtimes"))
  }

  # Sanity check sampling rates
  if (tr1@stats@sampling_rate != tr2@stats@sampling_rate) {
    stop(paste("correlationMetric: Incompatible sampling rates"))
  }
  
  # NOTE:  Correlation demands vectors of the same length (even though we will ignore NA).
  # NOTE:  We will truncate by up to one second to ensure this.
  
  # Sanity check lengths
  l1 <- length(tr1)
  l2 <- length(tr2)
  # Only complain if if sample lengths differ AND the difference is greater than the inverse sampling rate
  if ( (abs(l1 - l2) > 1) && (abs(l1 - l2) > 1/tr1@stats@sampling_rate) )  {
    stop(paste("correlationMetric: Incompatible lengths tr1 =",l1,", tr2 =",l2))      
  } else {
    min_length <- min(l1,l2)
  }

  # Sanity check network and station
  if (tr1@stats@network != tr2@stats@network || tr1@stats@station != tr2@stats@station) {
    stop(paste("correlationMetric: Incompatible trace ids '", tr1@id, "', '", tr2@id, "'", sep=""))
  }

  # Create two-channel ids if needed
  locations <- tr1@stats@location
  if (tr2@stats@location != tr1@stats@location) {
    locations <- paste(locations, ":", tr2@stats@location, sep="")
  }
  channels <- tr1@stats@channel
  if (tr2@stats@channel != tr1@stats@channel) {
    channels <- paste(channels, ":", tr2@stats@channel, sep="")
  }
  snclq <- paste(tr1@stats@network, tr1@stats@station, locations, channels, tr1@stats@quality, sep=".")

  # Calculate the correlation metric
  cor <- cor(tr1@data[1:min_length], tr2@data[1:min_length], use="na.or.complete")

  # Create and return a list of Metric objects
  m1 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="cross_talk", value=cor)

  return(c(m1))

}

