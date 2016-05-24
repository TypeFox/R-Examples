##
##    Metric to process a dataframe with daily means and return a vector of daily
##    likelihoods that a DC shift occurred. -- modified to return the daily likelihood
##	  for just the last day represented.
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
# REC - Jan 2014 - making modifications around the roll_sd call in order to address
# seg faults that are stopping production.
# REC - Nov 2014 - modifying code so that it only computes one value, representing the Nth day.
# TODO - refactor code so that it only operates on a single metric and doesn't calculate unnecessarily
# 		for an entire array
#
dailyDCOffsetMetric <- function(df, offsetDays=5,
                                outlierWindow=7,
                                outlierThreshold=6.0,
                                outlierSelectivity=0.1) {
  
  # Sanity checks
  if (class(df) != "data.frame") {
    stop(paste("dailyDCOffsetMetric: Argument 'df' is of class",class(df),". A dataframe is required."))
  }
  if (nrow(df) < 1) {
    stop(paste("dailyDCOffsetMetric: No data found in dataframe."))
  }
  if (!("sample_mean" %in% names(df))) {
    stop(paste("dailyDCOffsetMetric: Dataframe does not contain variable 'sample_mean'."))    
  }
  
  # Metric algorithm:
  #
  #   data0 = download length(df$sample_mean) of daily means
  #   data1 = remove outliers using MAD outlier detection with a (outlierWindow)-day window
  #   data2 = replace outliers with rolling median values using a (outlierWindow)-day window
  #   weights = calculate absolute lagged differences with 1-N day lags (big jumps have large values)
  #   metric0 = multiply the lagged differences together and take the N'th root
  #   stddev0 = calculate the rolling standard deviation of data2 with a N-day window
  #   METRIC = divide metric0 by the median value of stddev0

  # Replace outliers with rolling median values
  outliers <- seismicRoll::findOutliers(df$sample_mean,outlierWindow,outlierThreshold,outlierSelectivity,1)
  cleanMean <- df$sample_mean
  cleanMean[outliers] <- seismicRoll::roll_median(df$sample_mean,7)[outliers]
  
  # Vanilla metric
  metric <- rep(1.0,length(cleanMean))
  
  # Have a minimum value to prevent occasional zeros associated with different lags from completely wiping out large values
  diffMin <- rep(1e-3,length(cleanMean))

  # Create vectors of daily differences with N increasing lags, multiplying them together and then taking the N'th root
  for (i in seq(offsetDays)) {
    # Create a daily metric from the lagged data with NA's at the beginning.  Each date has the difference
    # between that date and the value 'i' days earlier.
    dailyDiff <- c(rep(NA,i),abs(diff(cleanMean,lag=i))) 
    # Multiplying them together weights those shifts that last for N days.
    metric <- metric * pmax(diffMin, dailyDiff) 
  }
  metric <- metric^(1/offsetDays)
  
  # Scale the metric by the median of the rolling sd with a window size of offsetDays
  #cat("DEBUG: seismicRoll::roll_sd cleanMean=", cleanMean, ", offsetDays=", offsetDays,"\n")
  #
  scaling <- stats::quantile(seismicRoll::roll_sd(cleanMean,offsetDays),0.5,na.rm=TRUE)
  metric <- metric / scaling
  
  # We have a vector of metric values and now convert these into a list of SingleValueMetric objects
  # REC -- we will bypass any looping and just take the last element
  metricList <- list()  # (kludge) this metricList will just have one element with the current implementation
  index <- length(metric) # index is fixed to the last element
  if (is.na(metric[index])) {
	  stop("dailyDCOffsetMetric: NA value found for dc_offset -- stop") 
  } else {
      metricList <- list ( new("SingleValueMetric",
                                 snclq=df[index,"snclq"], starttime=df[index,"starttime"], endtime=df[index,"endtime"],
                                 metricName="dc_offset", value=metric[index])   )
  }
  return(metricList)
}
