##
##    Metric containing state-of-health flag counts from the miniSEED record
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

stateOfHealthMetric <- function(st) {
  
  starttime <- st@requestedStarttime
  endtime <- st@requestedEndtime
  
  # Make sure we're working with a single snclq
  unique_ids <- uniqueIds(st)
  if (length(unique_ids) > 1) {
    stop(paste("basicStatsMetric: Stream has",length(unique_ids),"unique identifiers"))
  }
  snclq <- unique_ids[1]
  
  # act_flags
  # [1] Calibration signals present
  # [2] Time correction applied
  # [3] Beginning of an event, station trigger
  # [4] End of an event, station detrigger
  # [5] A positive leap second happened in this record
  # [6] A negative leap second happened in this record
  # [7] Event in progress
  # [8] Undefined bit set

  calibration_signal <- st@act_flags[1]
  timing_correction <- st@act_flags[2]
  event_begin <- st@act_flags[3]
  event_end <- st@act_flags[4]
  event_in_progress <- st@act_flags[7]
  
  # io_flags
  # [1] Station volume parity error possibly present
  # [2] Long record read (possibly no problem)
  # [3] Short record read (record padded)
  # [4] Start of time series
  # [5] End of time series
  # [6] Clock locked
  # [7] Undefined bit set
  # [8] Undefined bit set

  clock_locked <- st@io_flags[6]
                   
  # dq_flags
  # [1] Amplifier saturation detected
  # [2] Digitizer clipping detected
  # [3] Spikes detected
  # [4] Glitches detected
  # [5] Missing/padded data present
  # [6] Telemetry synchronization error
  # [7] A digital filter may be charging
  # [8] Time tag is questionable
  
  amplifier_saturation <- st@dq_flags[1]
  digitizer_clipping <- st@dq_flags[2]
  spikes <- st@dq_flags[3]
  glitches <- st@dq_flags[4]
  missing_padded_data <- st@dq_flags[5]
  telemetry_sync_error <- st@dq_flags[6]
  digital_filter_charging <- st@dq_flags[7]
  suspect_time_tag <- st@dq_flags[8]
  
  # timing_qual
  timing_quality <- st@timing_qual
  
  
  # Create and return a list of Metric objects
  m1 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="calibration_signal", value=calibration_signal)
  m2 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="timing_correction", value=timing_correction)
  m3 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="event_begin", value=event_begin)
  m4 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="event_end", value=event_end)
  m5 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="event_in_progress", value=event_in_progress)
                   
  m6 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="clock_locked", value=clock_locked)
                   
  m7 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="amplifier_saturation", value=amplifier_saturation)
  m8 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="digitizer_clipping", value=digitizer_clipping)
  m9 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="spikes", value=spikes)
  m10 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="glitches", value=glitches)
  m11 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="missing_padded_data", value=missing_padded_data)
  m12 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="telemetry_sync_error", value=telemetry_sync_error)
  m13 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="digital_filter_charging", value=digital_filter_charging)
  m14 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="suspect_time_tag", value=suspect_time_tag)
  m15 <- new("SingleValueMetric", snclq=snclq, starttime=starttime, endtime=endtime, metricName="timing_quality", value=timing_quality)
  
  return(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15))
  
}
