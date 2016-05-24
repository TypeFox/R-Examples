NULL
#' 
#' Calculates the average difference between the series of mean daily temperatures calculated by (Tmax + Tmin) / 2 and the average calculated by 24 hourly values a day, as resulting from the interpolation (or from measurements). The function works on data tables with series on columns.
#'
#' @title Calculates mean bias (difference between (max+min)/2 and 24-hour averages) in mean daily temperature series
#'
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param TMIN  data frame with daily minimum temperatures in columns. The first 3 columns are skipped (dates as year, month and day are supposed to be stored in these columns)
#' @param TMAX  same for TMAX
#' @param TMEAN same for TMEAN. Should come from 24-hour daily means.
#' @param min_valid  min nr. of valid days in a month for retaining its average value (if valid days are fewer, monthly value is \code{NA}). Default is 21.

#' @export 

#' @return A vector of means of daily biases, where the \code{TMEAN} is considered the "true" (reference) value

#' @references 
#' Eccel, E., 2010: What we can ask to hourly temperature recording. Part II: hourly interpolation of temperatures for climatology and modelling. Italian Journal of Agrometeorology XV(2):45-50 
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag45.pdf},\url{www.agrometeorologia.it}
#'
#' 
#'  See also: Eccel, E., 2010: What we can ask to hourly temperature recording. Part I: statistical vs. meteorological meaning of minimum temperature. Italian Journal of Agrometeorology XV(2):41-43.
#' \url{http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM\%202-2010_pag41.pdf},\url{www.agrometeorologia.it}
#' @note
#' Biases are calculated only on columns that are present in both \code{TMIN}/\code{TMAX} and \code{TMEAN}
#'

#' @examples
#' data(Trentino_hourly_T)
#' mo_bias <- bias(TMIN = Tn, TMAX = Tx, TMEAN = Tm_list, min_valid = 20)


#' @seealso \code{\link{daily_mean}}



###########################################################
# CHECKS THE BIAS BETWEEN 24-HOUR AND (max + min)/2 MEANS
###########################################################

bias<-function(TMIN, TMAX, TMEAN, min_valid=21)

{
# creates a data frame from the list
IDs_MEAN<-names(TMEAN[-1])
Tm_24.h<-TMEAN$Date
for(st in IDs_MEAN)
 Tm_24.h<-cbind(Tm_24.h, TMEAN[[st]])
names(Tm_24.h)[-(1:3)]<-IDs_MEAN

IDs_MIN.MAX<-names(TMIN[-(1:3)])
IDs<-setdiff(IDs_MEAN, setdiff(IDs_MEAN, IDs_MIN.MAX))
TEMPERATURE_AVERAGE<-(TMIN[,IDs] + TMAX[,IDs])/2
mm<-paste("0",Tm_24.h$month, sep=""); mm<-substr(mm, nchar(mm)-1, nchar(mm))
dd<-paste("0",Tm_24.h$day, sep=""); dd<-substr(dd, nchar(dd)-1, nchar(dd))
dates_MEAN<-paste(Tm_24.h$year, mm, dd)
mm_A<-paste("0",TMIN$month, sep=""); mm_A<-substr(mm_A, nchar(mm_A)-1, nchar(mm_A))
dd_A<-paste("0",TMIN$day, sep=""); dd_A<-substr(dd_A, nchar(dd_A)-1, nchar(dd_A))
TEMPERATURE_AVERAGE<-TEMPERATURE_AVERAGE[paste(TMIN$year, mm_A, dd_A) %in% dates_MEAN,]
bias<-data.frame(TEMPERATURE_AVERAGE - Tm_24.h[,IDs])
bias_annual<-colMeans(bias, na.rm=TRUE)
bias_annual[colSums(!is.na(bias)) < min_valid*12]<-NA
monthly_bias<-NULL
for(m in unique(mm))
{
  bias_mm<-bias[mm==m,]
  bias_means<-colMeans(bias_mm, na.rm=TRUE)
  bias_means[colSums(!is.na(bias_mm)) < min_valid]<-NA
  monthly_bias<-rbind(monthly_bias, bias_means)
}  
monthly_bias<-rbind(monthly_bias, bias_annual);row.names(monthly_bias)<-NULL
monthly_mean_bias <- rowMeans(monthly_bias, na.rm=TRUE);row.names(monthly_mean_bias)<-NULL
monthly_bias<-data.frame(monthly_bias, AVERAGE=monthly_mean_bias)
row.names(monthly_bias)<-c(unique(mm), "YEAR")
monthly_bias<-round(monthly_bias,2)

return(monthly_bias)

}
