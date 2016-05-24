NULL
#' 
#' Calibrates the monthly parameters for the definition of the interpolation algorithm values of: time of minimum temperature; time of maximum temperature; time of sunset; parameter 'c' (see references for details). It works on more than one series at once. It also calculates an average calibration table, for simulations of series with no calibration.
#'
#' @title Calibration of "hour" parameters for temperature interpolation
#' 
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param meas  measured hourly values file (table), where the first column is the series' ID
#' @param date.format date format in \code{meas}. Default is \code{"ymd"} (year month day)
#' @param cal_period calibration period for the series. Default is \code{NULL} (use whole period)
#' @param missing_value_code allows to deal missing data if a code (either numeric or character) is passed. Default is \code{NA}
#' @param min_valid_yrs minimum fraction of year(s) on any series ID necessary for carrying out calibration 
#' @param band_min (continuous) band of hours to seek day minimum 
#' @param band_max same for maximum time
#' @param band_suns same for sunset time
#' @param silent if set to \code{TRUE} removes notice of insufficient length for the calibration of single short series
#' @param aver_series set of series IDs (chr) used to calculate the average calibration (if \code{NULL}: all stations with valid data are included)
#' 
#' @export 

#' @return a list of calibration tables, one for each series, plus one average table at the bottom (named "Average")

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}
#'
#'  
#' Original algorithm from: Cesaraccio, C., Spano, D., Duce, P., Snyder, R.L., 2001. An improved model for determining degree-day values from daily temperature data. Int. J. Biometeorol. 45: 161-169.
#' \url{http://www.springerlink.com/content/qwctkmlq3tebthek/}
#'
#'  @note
#' \code{meas} must be organized as 4-field records, all series in the same file, no headers. Column order:  station ID, date, time (hour), T, [others fields, if any...] separated by spaces. This field order is mandatory.
#'
#' Default date format is "ymd" (yyyy/mm/dd). Different combinations can be passed to function with \code{date.format}, but separator must be "/"
#'
#' 
#' See also: Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.

#' @examples
#' data(Trentino_hourly_T)
#' stations <- c("T0001","T0010","T0129")
#' 
#' calibration_l <- par_calibration(meas = h_d_t[h_d_t$V1 %in% stations,], 
#'                  missing_value_code = -999.9, 
#' 					band_min = 4:8, band_max = 13:16, 
#' 					band_suns = 16:20, cal_period = NULL)

#' @seealso \code{\link{shape_calibration}}

############################################
# CALIBRATES PARAMETERS STATION BY STATION
############################################

par_calibration<- function(meas, date.format="ymd", cal_period=NULL, missing_value_code=NA, min_valid_yrs=1, band_min=4:9, band_max=12:16, band_suns=14:20, silent=FALSE, aver_series=NULL)

{

IDs_OK<- NULL
calibration_l<-NULL
IDs<-as.character(unique(meas[,1]))

for(sta in IDs)
{
if(silent==FALSE) 
  print(sta, quote=FALSE)
T_meas_hourly<-meas[meas[,1]==sta,] 
data<-date.mdy(as.date(as.character(T_meas_hourly[,2]), order=date.format))
T_meas_hourly$year<-data$year; T_meas_hourly$month<-data$month; T_meas_hourly$day<-data$day; T_meas_hourly$hour<-as.integer(substr(T_meas_hourly[,3], 1,2))
T_meas_hourly[,4][T_meas_hourly[,4]==missing_value_code]<-NA
T_meas_hourly<-data.frame(T_meas_hourly[,1], T_meas_hourly$year, T_meas_hourly$month, T_meas_hourly$day, T_meas_hourly$hour, T_meas_hourly[,4])
names(T_meas_hourly)<-c("ID", "year", "month", "day", "hour","T")

# filters the hourly series according to the calibration period
if(is.null(cal_period[1])) period<-c(min(T_meas_hourly$year, na.rm=TRUE), max(T_meas_hourly$year, na.rm=TRUE))  else {
  period<- c(cal_period[1], cal_period[length(cal_period)])
  if(cal_period[1] < min(T_meas_hourly$year)) {
    period[1]<-cal_period[1]; print("Warning: calibration start year forced to first year of hourly data", quote=FALSE)  }
  if(cal_period[length(cal_period)] > max(T_meas_hourly$year, na.rm=TRUE)) {
    period[2]<-cal_period[length(cal_period)]; print("Warning: calibration end year forced to last year of hourly data", quote=FALSE)  }
 } # else  
T_meas_hourly<-subset(T_meas_hourly, T_meas_hourly$year %in% period)

# calculates  the corresponding daily tables
T_band<-subset(T_meas_hourly, T_meas_hourly$hour>=min(band_min) & T_meas_hourly$hour<=max(band_min))
if(nrow(T_band) > sum(is.na(T_band$T)))
 Tn_period<-aggregate(T_band$T, by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=min) else 
 { year_month_day<-aggregate(T_band[,2:4], by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=unique, na.rm=TRUE)[4:6] 
  Tn_period<-data.frame(year_month_day,x=rep(NA, nrow(year_month_day))) }
Tmin<-data.frame(year=Tn_period$year, month=Tn_period$month, day=Tn_period$day, Tn=Tn_period$x)
T_band<-subset(T_meas_hourly, T_meas_hourly$hour>=min(band_max) & T_meas_hourly$hour<=max(band_max))
if(nrow(T_band) > sum(is.na(T_band$T)))
 Tx_period<-aggregate(T_band$T, by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=max) else 
 { year_month_day<-aggregate(T_band[,2:4], by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=unique, na.rm=TRUE)[4:6] 
   Tx_period<-data.frame(year_month_day,x=rep(NA, nrow(year_month_day))) }
Tmax<-data.frame(year=Tx_period[3], month=Tx_period[2], day=Tx_period[1], Tx=Tx_period$x)

Tsuns<-Tmin[1, 4] # Tsuns inizialized with Tn of the first day


# a check on the series length
if(nrow(T_meas_hourly[!is.na(T_meas_hourly$T),]) >=  365*24*min_valid_yrs)
{ # series length OK
IDs_OK<-c(IDs_OK, sta)


# calibrates Tmin

time_min_mode<-NULL  
matr_time_min<-NULL  

for(mm in 1:12)
{
h_count<-NULL
T_band<-subset(T_meas_hourly, T_meas_hourly$hour>=min(band_min) & T_meas_hourly$hour<=max(band_min) & T_meas_hourly$month==mm)
if(nrow(T_band) > sum(is.na(T_band$T)))
 Tn_band<-aggregate(T_band$T, by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=min) else
 { year_month_day<-aggregate(T_band[,2:4], by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=unique, na.rm=TRUE)[4:6] 
  Tn_band<-data.frame(year_month_day,x=rep(NA, nrow(year_month_day))) }
for(h in band_min)
 {
 T_h_band<-subset(T_band, T_band$hour==h)
 h_count[h-min(band_min)+1]<-nrow(subset(T_h_band,  T_h_band$T==Tn_band$x ))
 }
time_min_mode[mm]<-min(band_min[h_count==max(h_count)]) # in case of parity the earliest hour is chosen
matr_time_min<-cbind(matr_time_min, h_count)
}


# calibrates Tmax

time_max_mode<-NULL
matr_time_max<-NULL 	# solo per statistica (dalla 0 alle 23) 

for(mm in 1:12)
{
h_count<-NULL
T_band<-subset(T_meas_hourly, T_meas_hourly$hour>=min(band_max) & T_meas_hourly$hour<=max(band_max) & T_meas_hourly$month==mm)
if(nrow(T_band) > sum(is.na(T_band$T)))
 Tx_band<-aggregate(T_band$T, by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=max) else
{ year_month_day<-aggregate(T_band[,2:4], by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=unique, na.rm=TRUE)[4:6] 
  Tx_band<-data.frame(year_month_day,x=rep(NA, nrow(year_month_day))) }
for(h in band_max)
 {
 T_h_band<-subset(T_band, T_band$hour==h)
 h_count[h-min(band_max)+1]<-nrow(subset(T_h_band,  T_h_band$T==Tx_band$x ))
 }
time_max_mode[mm]<-max(band_max[h_count==max(h_count)]) # in case of parity the latest hour is chosen
matr_time_max<-cbind(matr_time_max, h_count)
}


# calibrates Tsuns (time of max T decrease)
# (NB: the result is decreased by 1 in the interpolation algorithm)
# creates a file of T differences (NB may generate "Warning" messages from R!)

delta.T<-NULL
delta.T[1]<-NA 
for(j in 2:nrow(T_meas_hourly))  {delta.T[j]<-T_meas_hourly$T[j]-T_meas_hourly$T[j-1]}
T_delta_h<-data.frame(T_meas_hourly[,1:5], delta=delta.T)

time_suns_mode<-NULL
for(mm in 1:12)
{
h_count<-NULL
T_band<-subset(T_delta_h, T_delta_h$hour>=min(band_suns) & T_delta_h$hour<=max(band_suns)  & T_delta_h$month==mm)
if(nrow(T_band) > sum(is.na(T_band$delta)))
 delta_band<-aggregate(T_band$delta, by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=min) else
 { year_month_day<-aggregate(T_band[,2:4], by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=unique, na.rm=TRUE)[4:6] 
  delta_band<-data.frame(year_month_day,x=rep(NA, nrow(year_month_day))) }
for(h in band_suns)
 {
 T_h_band<-subset(T_band, T_band$hour==h)
 h_count[h-min(band_suns)+1]<-nrow(subset(T_h_band,  T_h_band$delta==delta_band$x ))
 }
time_suns_mode[mm]<-min(band_suns[h_count==max(h_count)]) # in case of parity the latest hour is chosen
}


# calibrates "c"

# creates a daily file of T at sunset, knowing the average monthly sunset time
Tsuns_d<-subset(T_meas_hourly,T_meas_hourly$hour==time_suns_mode[T_meas_hourly$month])

# creates a daily series of Tn of the day after
Tn_after_d<-NULL

for(i in 1:(nrow(Tmin)-1))
 Tn_after_d[i]<-Tmin$Tn[i+1]

# the last day of all vectors except Tn_after_d must be removed
# check:
if( paste(Tmin$year,Tmin$month, Tmin$day)[1] != paste(Tsuns_d$year, Tsuns_d$month, Tsuns_d$day)[1] |
 paste(Tmin$year,Tmin$month, Tmin$day)[nrow(Tmin)] != paste(Tsuns_d$year, Tsuns_d$month, Tsuns_d$day)[nrow(Tsuns_d)]) 
 print(paste(sta, "Error in lengths!"), quote=FALSE) else  {
C_d<-data.frame(Tmin[-nrow(Tmin),1:3],c=(Tmax$Tx[-nrow(Tmax)] - Tsuns_d$T[-nrow(Tsuns_d)] ) / (Tmax$Tx[-nrow(Tmax)] - Tn_after_d))	
C_d$c[C_d$c==Inf | C_d$c==-Inf]<-NA }

if(nrow(C_d) > sum(is.na(C_d$c)))
 C_m<-as.data.frame(aggregate(C_d$c, by=list(month=Tmin$month[-nrow(Tmin)]),FUN=mean, na.rm=TRUE)) else
 { year_month_day<-aggregate(C_d[,2:4], by=list(day=C_d$day, month=C_d$month, year=C_d$year),FUN=unique, na.rm=TRUE)[4:6] 
  C_m<-data.frame(year_month_day,x=rep(NA, nrow(year_month_day))) }
C_m<-round(C_m$x,2) 


calibration<-data.frame(time_min=time_min_mode, time_max=time_max_mode, time_suns=time_suns_mode, C_m=C_m)
calibration_l<-append(calibration_l, list(calibration))

} # condition on series length 
else if(silent==FALSE) print(paste("Too missing data for ID", sta, "series"), quote=FALSE)



} # "sta in IDs"
names(calibration_l)<-IDs_OK

# creates a table with the average calibration


IDs_aver<-IDs_OK
if(!is.null(aver_series)) 
{
IDs_aver<-setdiff(aver_series, setdiff(aver_series,IDs_OK))
if(length(setdiff(aver_series, IDs_aver)) > 0)
 print(paste("Warning: series", setdiff(aver_series, IDs_aver), "excluded from average calibration table due too many missing data"), quote=FALSE)
}

time_min<-1:12  # month column
time_max<-1:12  
time_suns<-1:12 
C_m<-1:12
for(sta in IDs_aver) {
 time_min<-cbind(time_min, as.vector(as.data.frame(calibration_l[[sta]])[1]))
 time_max<-cbind(time_max, as.vector(as.data.frame(calibration_l[[sta]])[2]))
 time_suns<-cbind(time_suns, as.vector(as.data.frame(calibration_l[[sta]])[3]))
 C_m<-cbind(C_m, as.vector(as.data.frame(calibration_l[sta])[4]))  }
time_min<-time_min[-1]; names(time_min)<-IDs_aver
time_max<-time_max[-1]; names(time_max)<-IDs_aver
time_suns<-time_suns[-1]; names(time_suns)<-IDs_aver
C_m<-C_m[-1]; names(C_m)<-IDs_aver

calib_average<-data.frame(time_min=round(rowMeans(time_min), 0), time_max=round(rowMeans(time_max), 0), time_suns=round(rowMeans(time_suns), 0), C_m=round(rowMeans(C_m, na.rm=TRUE), 2))
calibration_l<-append(calibration_l, list(calib_average))

names(calibration_l)<-c(IDs_OK, "Average")

return(calibration_l)

}