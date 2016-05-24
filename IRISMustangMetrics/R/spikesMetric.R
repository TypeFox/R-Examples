##
##    Metric to calculate the number of outliers for a Stream of seismic data
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
# The IRIS DMC MUSTANG project keeps track of
#
#   - num_spikes
#
# This metric applies the roll_hampel filter from the 'seismicRolll' package.
# See:  "Median absolute deviation (MAD) outlier in Time Series"

spikesMetric <- function(st, windowSize=41, thresholdMin=6, selectivity=0.4) {
  
  starttime <- st@requestedStarttime
  endtime <- st@requestedEndtime
    
  # Calculate the total number of outliers from all traces
  count <- 0
  for (tr in st@traces) {
    
    result <- try( outlierIndices <- seismicRoll::findOutliers(tr@data, windowSize, thresholdMin, selectivity),
                   silent=TRUE )
    
    if (class(result)[1] == "try-error" ) {      
      err_msg <- geterrmessage()
      if (stringr::str_detect(err_msg,"n cannot be greater than length")) {
        # Skip very short trace segments
        next
      } else {
        ###MCRExitCode <- 1
        ###MCRWarning(paste("spikesMetric: skipping a trace in", tr@id, "with an unhandled error", err_msg))
        next
      } 
    } else {
      
      # NOTE:  Ignore adjacent outliers when determining the count of spikes.
      # NOTE:  But be sure there is at least one spike if there is at least one outlier.
      newCount <- length(which(diff(outlierIndices) > 1))
      if (length(outlierIndices) > 0 && newCount == 0) {
        newCount <- 1
      }
      count <- count + newCount
      
    }
    
  }  
  
  # Create and return a list of Metric objects
  snclq <- st@traces[[1]]@id
  m1 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="num_spikes", value=count)

  return(c(m1))
  
}
