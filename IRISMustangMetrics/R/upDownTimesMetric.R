##
##    Metric to extract the start- and end-times for all Traces in a Stream
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

# REC -- April 2014 - modified this routine to write out a different metric.  This one will write a simple metric for each
#   start and end time, representing channel uptime.

upDownTimesMetric <- function(st, min_signal=30, min_gap=60) {
  
  starttime <- st@requestedStarttime
  endtime <- st@requestedEndtime
  
  # Make sure we're working with a single snclq
  unique_ids <- uniqueIds(st)
  if (length(unique_ids) > 1) {
    stop(paste("meanMetric: Stream has",unique_ids,"unique identifiers"))
  }
  snclq <- unique_ids[1]

  # get the upDownTimes with a minimum gap specified in seconds
  # this is a vector of start,end,start,end....
  result <- try( upDownTimes <- getUpDownTimes(st, min_signal=min_signal, min_gap=min_gap),
                 silent=TRUE)

  # Handle error returns
  if (class(result)[1] == "try-error" ) {    
    # Write out the xml and stop with an error
    err_msg <- paste("ERROR in upDownTimesMetric(",snclq,",",min_signal,",",min_gap,"): ",geterrmessage(),sep="")
    stop(err_msg) 
  }    
  
  #OLD# Create and return a MultipleTimeValue metric
  #OLD#m1 <- new("MultipleTimeValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="up_down_times", values=upDownTimes)
  
  # NEW -- REC April 2014  
  # for each start/end pair, create separate metric entries
  up_time <- list()
  for (n in seq(1,length(upDownTimes)-1, by=2)) {
    # because of how n-1 samples are counted by the dataselect function
    # and how the dates are truncated to the second
    # we will take care of full-day cases by rounding the end time up to the next full hour
    # if both minutes and seconds are .gt. 59
    minute <- as.numeric(format(upDownTimes[n+1],"%M"))
    second <- as.numeric(format(upDownTimes[n+1],"%S"))
      if (minute >= 59 && second >= 59) {
        upDownTimes[n+1] <- round(upDownTimes[n+1],"mins")
      }
      
    # subtract start time from end time to get a duration value
    duration <- difftime(upDownTimes[n+1], upDownTimes[n], units="sec")
    if (is.na(duration)) {
      duration <- 0
    } else {
      intSeconds <- as.integer(duration)
    }
    m1 <- new("SingleValueMetric", snclq=snclq, starttime=upDownTimes[n], endtime=upDownTimes[n+1], metricName="channel_up_time", value=intSeconds)
    up_time <- append(up_time,m1)
  }

  return(c(up_time))
  
}
