#
# This file contains a script example to run a temperature interpolation with the package "Intepol.T"
#
# Author: Emanuele Eccel - 10-07-2012
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################


rm(list=ls())

library(Interpol.T)
data(Trentino_hourly_T)

#####################################################
# 1.
# CALIBRATION OF INTERPOLATION PARAMETERS
# HOURS OF MIN, MAX, SUNSET, COEFF. "c"
# (SKIP IF CALIBRATION TABLES ARE ALREADY AVAILABLE)
#####################################################

# set up parameters to be passed to function "par_calibration"
m_v_c<- -999.9   # missing value
c_p<-NULL        # calibration period
m_v_y<- 1        # min valid years for calibration
b_n<-4:9         # range of time for minimum
b_x<-12:16       # range of time for maximum
b_s<-14:20       # range of time for sunset
a_s<-c("T0147", "T0149", "T0150", "T0152", "T0154", "T0169","T0370") # series on which average calibration table is calculated


calibration_l<-par_calibration(meas=h_d_t, missing_value_code=m_v_c, min_valid_yrs=m_v_y, band_min=b_n, band_max=b_x, band_suns=b_s, date.format="ymd", silent=FALSE, aver_series=a_s, cal_period=c_p)


##############################################################
# 2. 
# CALIBRATION OF INTERPOLATION PARAMETERS
# DEFINITION OF ratio_dtr 
# (SKIP IF CALIBRATION TABLES ARE ALREADY AVAILABLE OR IF
# THE NIGHT-CURVE SHAPE CALIBR. IS NOT GOING TO BE USED)
##############################################################

ratio_dtr_r <- c(0,4)      # range of values for calibration of ratio_dtr

calibration_shape<-shape_calibration(meas=h_d_t, cal_times_list=calibration_l, band_min=0:23, band_max=0:23, ratio_dtr_range=ratio_dtr_r, min_mo.length=21, full.24.hrs.span_min=TRUE)

###########################################
# 3. 
# APPLICATION OF HOURLY INTERPOLATION
###########################################


start<-2004       # start of interpolation period
end<-2005         # end of interpolation period

Th_int_list<-Th_int_series(cal_times=calibration_l, TMIN=Tn, TMAX=Tx, start_year=start, end_year=end, full.24.hrs.span_min=TRUE)


###########################################
# 4. 
# GENERATION OF DAILY MEANS & CONTROL
###########################################


# generates the daily means from the hourly tables

Tm_list<-daily_mean(hourly_list=Th_int_list, series_names=NULL)


# calculates the bias between average of 24 hourly values and (Tmin + Tmax)/2

mo_bias<-bias(TMIN=Tn, TMAX=Tx, TMEAN=Tm_list, min_valid=20)

plot(mo_bias$AVERAGE[1:12], main= "(Tn+Tx)/2 - mean(24 values)", xlab="MONTH", ylab="DELTA T [degC]")
hist(mo_bias$AVERAGE[1:12], main= "Hist. of freq. (months) of delta.means", xlab="(Tn+Tx)/2 - mean(24 values)")

# plots charts for comparison between measured and simulated

m_v_c<- -999.9           # missing value code
start <- "1Jan2004"      # start date of charts
end <- "31Jan2004"       # end date of charts

plot_meas_sim(meas=h_d_t, sim= Th_int_list, missing_code=m_v_c, chart.start=start, chart.end=end, leg.pos="topright")
