##
##    Metric to calculate gaps and overlaps in a Stream of seismic data
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
#   - num_gaps
#   - max_gap
#   - num_overlaps
#   - max_overlap
#   - availability

gapsMetric <- function(st) {
  
  starttime <- st@requestedStarttime
  endtime <- st@requestedEndtime
  
  # Calcuate gaps
  gapInfo <- getGaps(st)
  
  # NOTE: getGaps.Stream() checks that all traces have the same ID
  
  if (sum(gapInfo$gaps) == 0) {
    num_gaps <- 0
    max_gap <- 0
    num_overlaps <- 0
    max_overlap <- 0
    gap_secs <- 0.0
  } else {
    # gaps
    gaps <- gapInfo$gaps[gapInfo$gaps > 0]
    num_gaps <- length(gaps)
    if (num_gaps > 0) {
      max_gap <- max(gaps, na.rm=TRUE)
    } else {
      max_gap <- 0
    }
    gap_secs <- sum(gaps[gaps>0])
    # overlaps
    overlaps <- gapInfo$gaps[gapInfo$gaps < 0]
    num_overlaps <- length(overlaps)
    if (num_overlaps > 0) {
      max_overlap <- abs(min(overlaps, na.rm=TRUE))
    } else {
      max_overlap <- 0
    }
    overlap_secs <- sum(gaps[gaps<0]) # overlap_secs is never used
  }
  
  if (num_gaps == 0) {
    percent_availability <- 100
  } else {
    totalSecs <- as.numeric(difftime(endtime,starttime,units="secs"))
    percent_availability <- 100 - 100 * gap_secs / totalSecs
  }
 
  # we do not want negative values, or values greater than 100, so contain this...
  if (percent_availability < 0) {
	  percent_availability <- 0
  }
  if (percent_availability > 100) {
	  percent_availability <- 100
  }
  
  snclq <- st@traces[[1]]@id

  # Create and return a list of Metric objects
  m1 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="num_gaps", value=num_gaps)
  m2 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="max_gap", value=max_gap)
  m3 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="num_overlaps", value=num_overlaps)
  m4 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="max_overlap", value=max_overlap)
  m5 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="percent_availability", value=percent_availability)

  return(c(m1,m2,m3,m4,m5))
  
}
