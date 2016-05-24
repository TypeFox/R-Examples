##
##    Metric to calculate the cross-correlation between two streams of seismic data
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
# Notes from Mary Templeton email:
#
#        Cross-correlate the current SNCL's data with the neighboring SNCL's data
#            Find the min and max of the cross-correlation function
#                Of these, save the value with the largest absoute value (preserve the sign)
#                Preserve the lag in seconds to the min or max value saved above
#                Find the difference in traveltime between the 2 traces: (current - neighbor)
#                Adjust the lag by subtracting this difference (preserving sign) from the lag
#
#        Report the 
#            current SNCL name (as usual)
#            neighbor SNCL name
#            largest absolute correlation function value (with origin sign) - we can call this peak_correlation
#            adjusted lag - we can call this peak_lag

### REC May 2014 -- we are now mapping the two measurements to the following metrics names
#
# peak_correlation --> polarity_check
# peak_lag --> timing_drift

crossCorrelationMetric <- function(st1, st2, maxLagSecs=10, filter=signal::butter(2,0.2)) {
  
  #-----------------------------------------------------------------------------
  # Sanity checks and data resampling
  #-----------------------------------------------------------------------------
  
  # Sanity check number of traces
  if (length(st1@traces) > 1) {
    stop(paste("crossCorrelationMetric:",st1@traces[[1]]@id,"has multiple traces."))
  }
  if (length(st2@traces) > 1) {
    stop(paste("crossCorrelationMetric:",st2@traces[[1]]@id,"has multiple traces."))
  }
  
  # Get the first (and only) trace and demean and detrend it 
  tr1 <- DDT(st1@traces[[1]],TRUE,TRUE,0)
  tr2 <- DDT(st2@traces[[1]],TRUE,TRUE,0)
  
  # Sanity check sampling rates
  if (tr1@stats@sampling_rate < 1) {
    stop(paste("crossCorrelationMetric:",tr1@id,"has a sampling_rate < 1."))
  }
  if (tr2@stats@sampling_rate < 1) {
    stop(paste("crossCorrelationMetric:",tr2@id,"has a sampling_rate < 1."))
  }
  
  # Deal with potentially different sampling rates
  sr1 <- as.integer(round(tr1@stats@sampling_rate))
  sr2 <- as.integer(round(tr2@stats@sampling_rate))
  sampling_rate <- min(sr1,sr2)
  
  if (sr1 > sampling_rate) {
    if (pracma::rem(sr1,sampling_rate) > 0) {
      stop(paste0("crossCorrelationMetric: sampling rates are not multiples of eachother:",
                  tr1@id,"=",sr1,", ",tr2@id,"=",sr2))
    }
    increment <- round(sr1/sampling_rate)
    d1 <- signal::decimate(tr1@data,increment)
  } else {
    d1 <- tr1@data
  }
  
  if (sr2 > sampling_rate) {
    if (pracma::rem(sr2,sampling_rate) > 0) {
      stop(paste0("crossCorrelationMetric: sampling rates are not multiples of eachother:",
                  tr1@id,"=",sr1,", ",tr2@id,"=",sr2))
    }
    increment <- round(sr2/sampling_rate)
    d2 <- signal::decimate(tr2@data,increment)
  } else {
    d2 <- tr2@data
  }
  
  # Sanity check that we have valid data everywhere
  if ( any(is.na(d1)) || any(is.na(d2)) ) {
    stop(paste("crossCorrelationMetric: NA values generated during resampling"))
  }
  
  #-----------------------------------------------------------------------------
  # Whew! Sanity checks are done. Now on to the metrics calculations.
  #-----------------------------------------------------------------------------
  
  # Apply low pass (or other) filter
  d1 <- signal::filter(filter,d1)
  d2 <- signal::filter(filter,d2)
  
  # Calculate cross-correlation
  lag.max <- sampling_rate * maxLagSecs 
  xcorr <- stats::ccf(d1, d2, lag.max, plot=FALSE)
  
  # Find strongest correlation 
  corrMin <- min(xcorr$acf)
  corrMax <- max(xcorr$acf)  
  if (abs(corrMin) > abs(corrMax)) {
    peak_correlation <- corrMin
  } else {
    peak_correlation <- corrMax
  }
  
  # Find the lag associated with the strongest correlation
  lagPoints <- xcorr$lag[ which(xcorr$acf == peak_correlation) ]
  peak_lag <- lagPoints * 1/sampling_rate
  
  # Create and return metrics
  m1 <- new("SingleValueMetric", snclq=tr1@id, starttime=tr1@stats@starttime, endtime=tr1@stats@endtime, 
            metricName="polarity_check", value=peak_correlation,
            attributeName="snclq2", attributeValueString=tr2@id)
  m2 <- new("SingleValueMetric", snclq=tr1@id, starttime=tr1@stats@starttime, endtime=tr1@stats@endtime,
            metricName="timing_drift", value=peak_lag,
            attributeName="snclq2", attributeValueString=tr2@id)
  
  return(c(m1,m2))
  
}

