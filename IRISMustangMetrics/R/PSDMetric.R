##
##    Metric to calculate the Power Spectral Density (PSD) for a Stream of 
##    seismic data.
##
##    Copyright (C) 2014  Mazama Science, Inc.
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
# This metric returns a list of average power spectra associated with a stream's
# data by breaking the Stream into chunks and calculating the spectrum for each.
#
# The 2004 McNamara paper "Ambient Noise Levels in the Continental United States"
# is available here:
# http://geohazards.cr.usgs.gov/staffweb/mcnamara/PDFweb/Noise_PDFs.html
#
#
# PSD algorithm from Mary Templeton on 2013-06-05
#   https://seiscode.iris.washington.edu/issues/46
#
# Target channels would include 
# {LMBSEH}{HNL}?
# {LB}G?
# 
# ----
# If channel is named L??
#   Z = 3 hours
# Else If channel is named M??
#   Z = 2 hours  (my recommendation, PQLX doesn't have a rule for M channels)
# Else
#   Z = 1 hour
# 
# Divide trace into Z-hour segments with 50% overlap
# Foreach Z-hour segment
#     Truncate segment to nearest power of 2 samples
#     Generate an averaged/smoothed PSD as follows
#         Divide each truncated Z-hour segment into 13 segments with 75% overlap (Z/4 seconds length each)
#         Foreach of 13 segments
#             Demean
#             Detrend
#             Apply 10% sine taper
#             FFT 
#             PSD 
#             Normalize the power at each PSD frequency, multiplying it by (2*dt/Nseg) where Nseg is the number of samples in the segment
#         Average the 13 resulting PSDs to get averaged PSD
#         Multiply averaged PSD by sine taper scale factor (= 1.142857)
#         Frequency smooth the averaged PSD over 1-octave intervals at 1/8-octave increments, 
#             (reducing # frequency samples by factor of 169; include 0.1 Hz as one of the geometric mean f values to sync f sampling)
#         Store the smoothed, averaged Z-hour PSD
# 


PSDMetric <- function(st,
                      expLoPeriod=4/(st@traces[[1]]@stats@sampling_rate), expHiPeriod=100,
                      linLoPeriod=16/(st@traces[[1]]@stats@sampling_rate), linHiPeriod=50) {
  
  # NOTE:  All the details about choosing window sizes, etc. are done in the psdList function
  
  # Use the psdList() function to apply the McNamara algorithm
  psdList <- psdList(st)
  
  # There is a remote possiblilty that different traces will have different quality identifiers
  snclqs <- sapply(psdList, getElement, "snclq")
  if (length(unique(snclqs)) > 1) {
    ids <- paste(unique(snclqs),collapse=",")
    stop(paste("PSDMetric: More than one SNCLQ in trace:",ids,sep=""))
  } else {
    snclq <- unique(snclqs)[1]
  }
  
  # Use PSDs to create SpectrumMetric objects
  spectrumMetricList <- list()
  index <- 1
  for (psd in psdList) {
    spectrumMetricList[[index]] <- new("SpectrumMetric", snclq=snclq, 
                                       starttime=psd$starttime, endtime=psd$endtime, metricName="psd",
                                       freqs=psd$freq, amps=psd$spec, phases=psd$freq*0)
    index <- index + 1
  }
  
  
  # Calculate SingleValueMetric to be stored in a separate list
  psdStats <- psdStatistics(psdList)
  
  # pct_above and pct_below are functions of frequency returned by psdStats
  avg_pct_above <- mean(psdStats$pct_above, na.rm=TRUE)
  avg_pct_below <- mean(psdStats$pct_below, na.rm=TRUE)
  
  # NOTE:  The dead_channel_exponential metric is calculated by fitting the PSD mean
  # NOTE:  line as seen in a PDF plot to an exponential and calculating the standard deviation of the residuals.
  # NOTE:  The mean of a healthy set of PSDs will have a very non-exponential shape and
  # NOTE:  large residuals while a "dead channel" will have a PSD mean that appears as
  # NOTE:  an exponential decay as a function of log10(period).
  # NOTE:
  # NOTE:  The dead_channel_exp metric looks for an overall exponential decay regardless
  # NOTE:  of the frequency range of the detector. Because of this, it lops off a few
  # NOTE:  frequency bands at either end regardless of the frequencies they represent.
  # NOTE:  As the frequency bands will be sensor specific because of the calculations in
  # NOTE:  the psdList() function, we can get by with just lopping an integer number of bands
  # NOTE:  rather than specifying the band pass region we want.
  
  # TODO:  Do we need to specify frequency bands for the dead channel metrics?
  
  # NOTE:  This algorithm is purely heuristic and resulted from a visual assessment
  # NOTE:  of PDF plots of channels known to be "dead".
  # NOTE:
  # NOTE:  Another metric fitting the PSD mean line as a linear function of log10(period)
  # NOTE:  over a specific band has been added. The parameters to this function 
  # NOTE:  allow users to modify the band over which fitting will occur.
  
  # Convert frequency to period
  period <- 1/psdStats$freq
  
  # Exponential fit metric:
  #   1) determine index range (inverted because of freq->period conversion)
  #   2) Convert band from PSD mean to positive, non-zero values
  #   3) Fit log10(PSD mean) vs. log10(period) to a line
  #   4) Calculate the standard deviation of residuals of the fit (How close to exponential is the PSD mean line?)
  first <- max(which(period >= expHiPeriod)) + 1
  last <- min(which(period <= expLoPeriod)) - 1
  positiveMean <- psdStats$mean[first:last] - min(psdStats$mean[first:last]) + .1  
  expFit <- stats::lm(log10(positiveMean) ~ log10(period[first:last]))
  dead_channel_exp <- sd(expFit$residuals)
  
  # Create metricList
  starttime <- st@traces[[1]]@stats@starttime
  endtime <- st@traces[[length(st@traces)]]@stats@endtime
  
  m1 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="pct_above_nhnm", value=avg_pct_above)
  m2 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="pct_below_nlnm", value=avg_pct_below)
  m3 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="dead_channel_exp", value=dead_channel_exp)
  
  # Put these new metrics in svmList
  svMetricList <- list(m1,m2,m3)
  
  # NOTE:  Here we do something atypical.
  # NOTE:  Instead of returning a metricList, we return a list of lists.  
  listOfLists <- list(svMetricList=svMetricList,
                      spectrumMetricList=spectrumMetricList)
  
  
  return(listOfLists)
  
}

