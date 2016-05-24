NULL
#' @description Creates climate mean monthly values from a monthly series of temperature and precipitation.
#' 
#' @param series the monthly series of temperature and precipitation. 
#' @param first.yr first year of the period over which climatology is calculated
#' @param last.yr last year of the period over which climatology is calculated
#' @param max.perc.missing maximum acceptable percentage of missing data in the averaging period from \code{first.yr} to \code{last.yr} (0-99).
#'
#' @title Climate normals
#' @author Emanuele Eccel
#'
#' @return A data frame with climatic monthly values of: precipitation, minimum and maximum temperatures (if existing in series), mean temperature (either averaged from existing values in series, or calculated by the function as (Tn + Tx)/2), absolute minimum monthly temperature.
#'
#' @details 
#' \code{series} is a data frame with years, months, temperature (and precipitation) values. Names in series columns must include: year, month, Tn and Tx (minimum and maximum temperatures, respectively) or, as an alternative, Tm (mean temperatures).
#' 
#' 
#' If \code{first.yr} or  \code{last.yr} are NULL (default), the lowest and highest values in series are taken as the period.
#' 
#' @importFrom stats aggregate
#' @export
#' 
#' 
#' @examples
#' 
#' data(Trent_climate)
#' 
#' # clima_81_10 is a list of data frames of the type series, 
#' # each one referring to one station 
#' # having climatic means of temperature and precipitation 
#'
#' clima_81_10<-lapply(lista_cli, FUN=climate, first.yr=1981, last.yr=2010, max.perc.missing=15)
#'
climate <- function(series, first.yr=NULL, last.yr=NULL, max.perc.missing)
{
  if(is.null(first.yr)) first.yr <- min(series$year)
  if(is.null(last.yr)) last.yr <- max(series$year) 
  series_period<-series[series$year>=first.yr & series$year<=last.yr,]  
  series_cli_med<-aggregate(series_period, by=list(series_period$month), FUN=mean, na.rm=T)[-(1:2)]
  if(sum(!is.na(series_period$Tn)) >0 & "Tn" %in% names(series_period))
    series_abs_Tn<-aggregate(data.frame(series_period$month, series_period$Tn), by=list(series_period$month), FUN=min, na.rm=T)[3] else
    series_abs_Tn<-as.numeric(rep(NA,12))
  names(series_abs_Tn)<-"AbsTn"
  missing<-aggregate(series_period, by=list(series_period$month), FUN=function(x) 
    {count<-sum(is.na(x))
     return(count)} )[-(1:3)]
  series_cli_med[-1][missing > max.perc.missing/100*(last.yr - first.yr +1)] <- NA
  if("Tn" %in% names(series) & "Tx" %in% names(series) & !"Tm" %in% names(series)) 
    series_cli<-round(data.frame(series_cli_med, Tm=(series_cli_med$Tn + series_cli_med$Tx)/2, AbsTn=series_abs_Tn), 1) else
      series_cli<-round(data.frame(series_cli_med, AbsTn=series_abs_Tn), 1)
  if("Tn" %in% names(series)) series_cli$AbsTn[is.na(series_cli$Tn)]<-NA
  return(series_cli)
}
