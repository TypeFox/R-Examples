NULL
#' 
#' Calibrates the shape of the night interpolating curve, either horizontal-axe parabola or line, by changing the exponent z (see reference). It functions according to the comparison of the daily thermal range and the climate (reference) monthly one. 
#'
#' @title Calibrates the shape of the night interpolating curve
#' 
#' 
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param meas  measured hourly values file (table), where the first column is the series' ID
#' @param date.format input date format (formats for function \code{chron})
#' @param cal_times_list calibration list of "time" parameters (output of \code{\link{par_calibration}})
#' @param band_min band of hours of occurrence of day minimum in the daily series (continuous). See \code{Note}
#' @param band_max same for maximum time
#' @param ratio_dtr_range range for seeking the optimal value of \code{ratio_dtr}
#' @param nr_cycles number of calibration trials within the calibration ranges (all)
#' @param min_mo.length minimum number of days to calculate any monthly values of dtr (is passed to function \code{\link{Mo.Th.Ra.}})
#' @param full.24.hrs.span_min logical, if set to \code{FALSE} does not allow to shift minimum time to the late hours of the day
#' @param silent logical, if set to \code{TRUE} suppresses any warning issue
#' 

#' @export 

#' @return A list containing the optimum values of \code{ratio_dtr} 

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}
#' 
#' 
#' Original algorithm from: Cesaraccio, C., Spano, D., Duce, P., Snyder, R.L., 2001. An improved model for determining degree-day values from daily temperature data. Int. J. Biometeorol. 45: 161-169.
#' \url{http://www.springerlink.com/content/qwctkmlq3tebthek/}
#'  
#' 
#' See also: Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf},\url{www.agrometeorologia.it}
#' 
#' @note
#' \code{meas} must be organized as 4-field records, all series in the same file, no headers. Column order:  station ID, date, time (hour), T, [others fields, if any...] separated by spaces. This field order is mandatory.
#'
#' Default date format is "ymd" (yyyy/mm/dd). Different combinations can be passed to function with \code{date.format}, but separator must be "/"
#'
#' \code{band_min} and \code{band_max} are the time bands according to which the minimum and maximum temperature were calculated in the daily series to be interpolated. In general, they range from 0 to 23, unless the series has had some restriction in the calculation of minimum and maximum values. Hence, these bands can be different from those used to calibrate the most frequent occurrence of min and max.
#'
#' The optimal value of \code{ratio_dtr} (\code{k}, eq. 7, in the quoted reference Eccel (2010a)) is chosen as the one with the (absolute) minimum value of the bias (irrespective of its sign). \code{ratio_dtr} is the ratio between the daily thermal range of the day to be interpolated and the mean monthly value for that series. The corresponding values of mean absolute error and RMSE can be checked in the resulting list. 
#'
#' \code{min_mo.length} (passed to function \code{Mo.Th.Ra.}) refers to the sum of days along all the series for each specific month, not for any single month in one year.
#'
#' \code{full.24.hrs.span_min} is \code{TRUE} as default. If must be set to \code{FALSE} only if minimum values of the daily series have been calculated on a restricted time band, which is included in \code{band_min} (see function \code{par_calibration}).
#' If this is the case, the minimum of the interpolated curve will always fall within \code{band_min} (early hours of the day). If this option is erroneously chosen, errors as large as 0.6  deg C can arise in the average of mean daily T.


#' @examples
#' library(Interpol.T)
#' data(Trentino_hourly_T)
#' 
#' stations <- c("T0001","T0010","T0129")
#' 
#' calibration_shape <- shape_calibration(meas = h_d_t[h_d_t$V1 %in% stations,], 
#' 						cal_times_list = calibration_l[stations], 
#'                      band_min = 0:23, band_max = 0:23, ratio_dtr_range = c(0,4), 
#' 						min_mo.length=21)


#' @seealso \code{\link{par_calibration}}, \code{\link{Th_interp}}


###################################################
# CALIBRATES ratio_dtr 
###################################################


shape_calibration<- function(meas, date.format="ymd", cal_times_list, band_min=0:23, band_max=0:23, ratio_dtr_range=c(0,6),  nr_cycles=10, min_mo.length=21, full.24.hrs.span_min=TRUE,  silent=FALSE)

{

options(warn=-1)             # suppresses R warnings

IDs_OK<-names(cal_times_list)                     # all except "Average"
IDs_OK<- IDs_OK[IDs_OK!="Average"]
T_h_all<-meas[,1:4]
names(T_h_all)<-c("ID","date", "hour","T")
IDs<-as.character(unique(T_h_all$ID))                            # all stations IDs

delta_k.ratio<- (ratio_dtr_range[2]- ratio_dtr_range[1])/(nr_cycles-1)
ratio_dtr<- ratio_dtr_range[1] + (0:(nr_cycles-1))*delta_k.ratio

calibration_shape_ratio<-NULL

if(silent==FALSE)  {
print("Calibration of ratio_dtr...", quote = FALSE)
print("", quote=FALSE)  }


for(sta in IDs)
{

if(!silent) print(sta, quote=FALSE)

E.ratio<-NULL
E.abs.ratio<-NULL
RMSE.ratio<-NULL
E.lim<-NULL
E.abs.lim<-NULL
RMSE.lim<-NULL

# if sta OK loads the corresponding calibr. table and substit. missing parameters with the average ones; if sta not OK chooses the average table
if(sta %in% IDs_OK) {
calibr<-cal_times_list[[sta]] 
calibr[is.na(calibr)]<- cal_times_list$Average[is.na(calibr)] } else calibr<-cal_times_list$Average

T_meas_hourly<-T_h_all[T_h_all$ID==sta,]

data<-date.mdy(as.date(as.character(T_meas_hourly$date), order=date.format))
T_meas_hourly<-data.frame(T_meas_hourly$ID, data$year, data$month, data$day, T_meas_hourly$hour, T_meas_hourly$T)
names(T_meas_hourly)<-c("ID", "year", "month", "day", "hour","T")

# calculates the corresponding daily tables
# bands are modified because both min and max are generally sought over the whole 24-hour domain
T_band<-subset(T_meas_hourly, as.integer(T_meas_hourly$hour)>=min(band_min) & as.integer(T_meas_hourly$hour)<=max(band_min))
Tn_period<-aggregate(T_band$T, by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=min,na.rm=TRUE)
Tn_period$x[Tn_period$x==Inf]<-NA
Tmin<-data.frame(year=Tn_period[3], month=Tn_period[2], day=Tn_period[1], T=Tn_period$x)
T_band<-subset(T_meas_hourly, as.integer(T_meas_hourly$hour)>=min(band_max) & as.integer(T_meas_hourly$hour)<=max(band_max))
Tx_period<-aggregate(T_band$T, by=list(day=T_band$day, month=T_band$month, year=T_band$year),FUN=max,na.rm=TRUE)
Tx_period$x[Tx_period$x==-Inf]<-NA
Tmax<-data.frame(year=Tx_period[3], month=Tx_period[2], day=Tx_period[1], T=Tx_period$x)

Tsuns<-Tmin[1, 4] # initializes by approximating first Tsuns with the min of the same day

# calculates monthly means of daily thermal range
dtr<-Mo.Th.Ra.(Tmin=Tmin, Tmax=Tmax, name=sta, min_mo.length=min_mo.length, silent=silent)  # function for calculation of monthly thermal range

for(k in 1:nr_cycles)
{

Th_simul<-NULL  
for(day in 1:(nrow(Tmin)-1)) 
 {  
 if(exists("Th_l")) 
 { Ts<-Th_l$Tsuns; T24_bf <- Th_l$Th_24_before } else {Ts<-NULL; T24_bf<-NULL}
 Th_l<-Th_interp(Tmin=Tmin, Tmax=Tmax, Tsuns=Ts, Th_24_before=T24_bf, day=day, tab_calibr=calibr, dtr_month=dtr, ratio_dtr=ratio_dtr[k], late_min=full.24.hrs.span_min)
 Th<-Th_l$Th
 Th_simul<-append(Th_simul, round(Th,1))
 }

comparison<-data.frame(T_meas_hourly[1:length(Th_simul),], simul=round(Th_simul,1)); names(comparison)<-c("station", "year", "month", "day", "hour", "T.meas", "T.sim")
comparison<-data.frame(comparison, error=comparison$T.sim - comparison$T.meas, abs.err=abs(comparison$T.sim - comparison$T.meas))

e<-mean(comparison$error,na.rm=TRUE)
e.abs<-mean(comparison$abs.err,na.rm=TRUE)
rmse<-(sum(comparison$error^2,na.rm=TRUE)/sum(!is.na(comparison$error)))^0.5

E.ratio[k]<-e; E.abs.ratio[k]<-e.abs; RMSE.ratio[k]<-rmse

}   #  k cycle

# finds ratio_dtr that yields min errors and adds abs.errors and rmse 
# and creates error tables for all stations in IDs

min_errors_ratio<-data.frame(ratio.min=NA, ratio_err.min=NA, ratio_abs.min=NA, ratio_rmse.min=NA)
if(sum(is.na(E.ratio))!=nr_cycles)  # enough data for calibration
{
k.min<-which(abs(E.ratio)==min(abs(E.ratio)))[1]  # chooses (arbitrarily) the first, if more than 1
ratio.min<-ratio_dtr[k.min]
ratio_err.min<-E.ratio[k.min] 
ratio_abs.min<-E.abs.ratio[k.min]
ratio_rmse.min<-RMSE.ratio[k.min]
min_errors_ratio<-round(data.frame(ratio.min=ratio.min, ratio_err.min=ratio_err.min, ratio_abs.min=ratio_abs.min, ratio_rmse.min=ratio_rmse.min),3)
}
row.names(min_errors_ratio)<-sta

calibration_shape_ratio<-rbind(calibration_shape_ratio, min_errors_ratio)

} # stations

calibration_shape<-list(ratio=calibration_shape_ratio)

if(!silent) {
print("",quote=FALSE)
print("Calibration of ratio_dtr completed successfully!", quote=FALSE)
print("",quote=FALSE) }

calibration_shape<-list(ratio=calibration_shape_ratio)

options(warn=0)

return(calibration_shape)

}

