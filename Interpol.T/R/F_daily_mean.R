NULL
#' 
#' The function works on a list of hourly temperature series. The hourly list is the output of the interpolation function \code{Th_int}, called iteratively to produce a list where each component represents one interpolated series.
#'
#' @title Production of daily means from hourly temperature series
#'
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param hourly_list  the list of hourly temperatures
#' @param series_names  names of the serie to be averaged (if \code{NULL} (default): all series)

#' @export 

#' @return A list of daily averaged series

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#'
#' 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.
#' 
#' 
#' @note
#' The first element of \code{hourly_list} must be a data frame named "Date" and its columns "year", "month", "day" (a fourth column ("hours") is not used in this function)

#' @examples
#' data(Trentino_hourly_T)
#' # generates daily means for series T0001 and T0129:
#' Tm_list <- daily_mean(hourly_list = Th_int_list, series_names = c("T0001", "T0129"))

#' @seealso \code{\link{Th_int_list}}


###############################################################
# CREATES DAILY MEAN SERIES OF TEMPERATURE FROM HOURLY SERIES
###############################################################



daily_mean<-function(hourly_list, series_names=NULL)

{
daily_list<-NULL
IDs<-names(hourly_list)[-1]
if(!is.null(series_names))  IDs<-series_names   
 if(sum(series_names %in% names(hourly_list)[-1]) != length(series_names) ) print("Error: input name(s) not matching hourly series!", quote=FALSE) else
{

 date_g<-aggregate(hourly_list$Date, by=list(hourly_list$Date$day, hourly_list$Date$month, hourly_list$Date$year), FUN = unique)[,4:6]
 daily_list<-list(date_g)

 for(st in IDs)
 {
 print(paste("Now processing", st), quote=FALSE)
 serie_g<-data.frame(Tm=round(aggregate(hourly_list[[st]], by=list(hourly_list$Date$day, hourly_list$Date$month, hourly_list$Date$year), FUN=mean)[,4], 1))
 daily_list<-append(daily_list, list(serie_g))
 }
 names(daily_list)<-c("Date", IDs)
 }

return(daily_list)

}