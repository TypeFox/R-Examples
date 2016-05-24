NULL
#' 
#' Creates a table of mean monthly DTR (difference between maximum and minimum temperature) over the series' period, to use a reference for establishing the exponent "z" of the night stretch of the interpolation curve.
#'
#' @title Creates a table of mean monthly "daily thermal range" (dtr)
#'
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param Tmin  a table of 4 named columns, the first are year, month, day, the 4th is minimum temperature. The column names "month" and "T" are mandatory
#' @param Tmax  same for Tmax
#' @param name  name of the series dtr is calculated for
#' @param min_mo.length  minimum number of days necessary to calculate monthly dtr values
#' @param silent  logical. If \code{TRUE} no warning is issued

#' @export 

#' @return A vector of 12 dtr values (1 for January, ... 12 for December)


#' @examples
#' data(Trentino_hourly_T)
#'id<-"T0001"
#'Tmin<-data.frame(Tn[,1:3], T=Tn[,id])
#'Tmax<-data.frame(Tx[,1:3], T=Tx[,id])
#' dtr <- Mo.Th.Ra.(Tmin = Tmin, Tmax = Tmax, name = id)


#' @seealso \code{\link{shape_calibration}}, \code{\link{Th_int_series}}



######################################
# CALCULATES MONTHLY THERMAL RANGES
# CHECKING MISSING DATA
######################################

Mo.Th.Ra.<-function(Tmin, Tmax, name, min_mo.length= 21, silent=FALSE)

{

mens_Tn<-aggregate(Tmin, by=list(Tmin$month), FUN=mean, na.rm=TRUE); mens_Tx<-aggregate(Tmax, by=list(Tmax$month), FUN=mean, na.rm=TRUE)
dtr<-(mens_Tx - mens_Tn)[,5]
# cuts dtr values calculated with too few data
Tmin_nna<-Tmin[!is.na(Tmin$T),]
if(nrow(Tmin_nna) > 0) {
 mo.length<-aggregate(Tmin_nna$T, by=list(Tmin_nna$month, Tmin_nna$year), FUN=length)
 mo.length.mean<-aggregate(mo.length$x, by=list(mo.length$Group.1), FUN=sum, na.rm=TRUE)
 dtr[mo.length.mean$x < min_mo.length]<-NA  }

# checks missing terms and substitutes with previous months where possible
dtr_corr<-dtr
for(j in 2:12)
 if(is.na(dtr[j]) & !is.na(dtr[j-1])) {dtr_corr[j]<-dtr_corr[j-1]; if(silent==FALSE) print(paste(name, "Warning: dtr for month", j, "forced to dtr for month", j-1, "due to lack of valid data"), quote=FALSE)}
if(is.na(dtr[1]) & !is.na(dtr[12])) {dtr_corr[1]<-dtr_corr[12]; if(silent==FALSE) print(paste(name, "Warning: dtr for month 1 forced to dtr for month 12 due to lack of valid data"), quote=FALSE)}
dtr<-dtr_corr
if(sum(is.na(dtr))!=0)  # still missing terms: tries with successive months
for(j in 1:11)
 if(is.na(dtr[j]) & !is.na(dtr[j+1])) {dtr_corr[j]<-dtr_corr[j+1]; if(silent==FALSE) print(paste(name, "Warning: dtr for month", j, "forced to dtr for month", j+1, "due to lack of valid data"), quote=FALSE)} 
dtr<-dtr_corr
if(sum(is.na(dtr))!=0)  # still missing terms: no other substitution possible
 if(silent==FALSE) print(paste(name, "Warning: too few valid data for estimation of missing dtr terms"), quote=FALSE)

return(dtr)

}