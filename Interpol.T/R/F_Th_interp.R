NULL
#' 
#' The function creates 24 values of hourly temperature from minimum and maximum daily values. This function applies to single series and to single day couples of minimum and maximum temperature. It is called by functions \code{Th_int_series} and \code{shape_calibration}.
#' The function uses four different curves: from time 00 to the minimum time: a horizontal-axis parabola (a line, if this choice is enabled and according to the daily thermal range of the day); from minimum to maximum time: an increasing sinusoidal curve; from maximum time to sunset: a decreasing sinusoidal curve; 
#' from sunset to time = 23: a horizontal-axis parabola (a line, if this choice is enabled and according to the daily thermal range of the day).
#' Calibration parameters are series- and monthly-specific.
#' This function is operationally called by \code{Th_int_series}, which requires the daily series and the calibration table as input (plus other parameters). A general user will conveniently use the latter function. 
#'
#' @title 24-hourly interpolation of temperature
#' 
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param Tmin  a daily table of 4 named columns, the first 3 being year, month, day, the 4th minimum temperature. The column names "month" and "T" are mandatory
#' @param Tmax  same for Tmax
#' @param Tsuns  temperature at sunset time
#' @param Th_24_before  temperature at time 24 of the previous day (time 00 of the present day)
#' @param day  progressive number of the day (row of both \code{Tmin} and \code{Tmax}), corresponding to a day
#' @param tab_calibr  "hour" parameter calibration table for the specific series. See \code{\link{par_calibration}}
#' @param dtr_month  monthly daily thermal range table (see function \code{Mo.Th.Ra.})
#' @param ratio_dtr  parameter for the choice of the night curve shape; it is \code{NULL} if no calibration_shape is passed to the function by \code{Th_int_series}
#' @param late_min  logical; allows to shift the time of occurrence of minima to the late hours of the day (assumes the value of \code{full.24.hrs.span_min}, passed by functions \code{Th_int_series)} and \code{shape_calibration})

#' @export 

#' @return A vector containing the values from hour = 0 (element 1) to hour = 23 (element 24)

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}
#' 
#' 
#' Original algorithm from: Cesaraccio, C., Spano, D., Duce, P., Snyder, R.L., 2001. An improved model for determining degree-day values from daily temperature data. Int. J. Biometeorol. 45: 161-169.
#' \url{http://www.springerlink.com/content/qwctkmlq3tebthek/}
#'  	
#' See also: Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf},\url{www.agrometeorologia.it}
#' 
#' 
#' @note
#' The function is called by \code{Th_int_series}.
#'
#' If the series ID coincides with one with non-null results of the \code{par_calibration} function (enough data for calibration) its table is passed to the interpolation function, otherwise the average (\code{cal_table}) is used.
#'
#' A non-NULL value for \code{ratio_dtr} enables the function to interpolate night values with a line, if the conditions on Daily Thermal Range occur. This may give rise to a sharp change for the hours following the min.
#'
#' Tmin of the day before the first is set = to Tmin of the first day and Tmin of the day after the last = Tmin of the last day.
#'
#' If T from sunset falls below the minimum of the day (temptatively attributed to \code{time_min}), an adjustement is done: the early hours assume a constant T = T[00] and the minimum is shifted so that T[23] = \code{Tmin} for that day
#'
#' Since the very first value of T series at sunset (of the day before) is \code{NULL}, the first hourly values produced till \code{time_min} are =  \code{Tmin} of the day.


#' @seealso \code{\link{Th_int_list}}



#########################################################
# INTERPOLATES DAILY TEMPERATURE (ONE DAY, ONE SERIES)
# CALLED BY Th_int_series AND BY shape_calibration
#########################################################


Th_interp<-function(Tmin, Tmax, Tsuns=NULL, Th_24_before=NULL, day, tab_calibr, dtr_month=NULL, ratio_dtr=NULL, late_min=TRUE)

{

options(warn=-1)
Th<-NULL
mm<-Tmin$month[day]
Tn<-Tmin$T[day] 
Tx<-Tmax$T[day] 
if(is.null(Tsuns))
 Tsuns_d.before <- Tn + 1.0 else Tsuns_d.before<-Tsuns
if(is.null(Th_24_before))
 Th_24_before <- Tn
if(!is.na(Tmin[day+1,4])) Tn_after<-Tmin$T[day+1] else Tn_after<-Tmin$T[day]
if(!is.null(ratio_dtr)) {
 cloudy_morning<-!is.na(Tx-Tn) & (Tx-Tn)/dtr_month[Tmin$month[day]]<=ratio_dtr  
 cloudy_night<-!is.na(Tx-Tn_after) & (Tx-Tn_after)/dtr_month[Tmin$month[day]]<=ratio_dtr  }
if(late_min==FALSE) # Tmin has been calculated (and must fall) at early hours
 Tsuns<- Tx - tab_calibr$C_m[mm]*(Tx-Tn_after) else  # Tmin cal also fall at late hours
 Tsuns<- min(max(Tx - tab_calibr$C_m[mm]*(Tx-Tn_after), Th_24_before+2.0, na.rm=TRUE), Tx, na.rm=TRUE)
if(!is.na(Tsuns) & !is.null(Tsuns) & (Tsuns==Inf | Tsuns== - Inf)) Tsuns<-NA
time_min<-tab_calibr$time_min[mm]
time_max<-tab_calibr$time_max[mm]
time_suns<-tab_calibr$time_suns[mm]

# begins interpolation

# h = 00 to min
z<-0.5 
Tn_original<-Tn
if(!is.null(ratio_dtr))
 if(cloudy_morning==TRUE) z<-1 
b<-(Tn-Th_24_before)/(time_min^z)
for(h in 0:(time_min))  #  Th[1] = Th at time h = 0 and so on
 Th[h+1]<-Th_24_before + b*((h)^z)

# h = min to max
for(h in time_min:time_max) 
 Th[h+1]<-Tn+(Tx-Tn)/2 * ( 1 + sin((h-time_min)/(time_max-time_min)*pi - pi/2)) 

# h = max to sunset
for(h in (time_max):(time_suns)) 
 Th[h+1]<- Tsuns + (Tx-Tsuns)*sin(pi/2*(1+(h-time_max)/(time_suns-time_max)))     

# h = sunset to 23
z<-0.5 
if(!is.null(ratio_dtr))
 if(cloudy_night==TRUE) z<-1 
b<-(Tn_after - Tsuns) / ((time_min + 24-time_suns)^z)
for(h in time_suns:23)
 Th[h+1]<-Tsuns + b*((h-time_suns)^z)
Th_24_before<-Tsuns + b*((24-time_suns)^z)

# checks if T[23] is lower than Tmin: in this case, re-calculates the last period
if(late_min == TRUE & !is.na(Th[24] < Tn) &  Th[24] < Tn)
{
b<-(Tn - Tsuns) / ((23-time_suns)^z) 
for(h in time_suns:23)
 Th[h+1]<-Tsuns + b*((h-time_suns)^z)
Th_24_before<-Tsuns + b*((24-time_suns)^z)

# h = 00 to min, new
Th[1:time_min]<-Tsuns_d.before

# h = min to max, new
for(h in time_min:time_max) 
 Th[h+1]<-Th[time_min]+(Tx-Th[time_min])/2 * ( 1 + sin((h-time_min)/(time_max-time_min)*pi - pi/2)) 

}

Th_list<-list(Th,Tsuns,Th_24_before); names(Th_list)<-c("Th", "Tsuns","Th_24_before")
return(Th_list)

}
