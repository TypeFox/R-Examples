NULL
#' 
#' Generates dates and hours between start and end terms. The result is a data frame with years, months, days and hours. Called in function \code{\link{Th_int_series}}
#'
#' @title Generation of dates and hours between start and end terms 
#' 
#' 
#' @author  Emanuele Eccel, Emanuele Cordano \email{emanuele.eccel@@iasma.it}
#' 
#' @param day.begin begin date (day) - format specified in \code{date.format}
#' @param day.end end date (day) - format specified in \code{date.format}
#' @param h.begin begin time (hour - integer)
#' @param h.end end time (hour - integer)
#' @param date.format input date format (formats for function \code{chron})
#' 


#' @export 

#' @return A 4-column table having the folloging fields (integer): "year", "month", "day", "hour"

#' @note
#' Input start and end dates as \code{character} (default format "yyyy/mm/dd"), hours as \code{integers} (0 to 23).
#'
#' Date format can be changed according to package \code{chron}'s standard, e.g.  "y/m/d" (default) or "m/d/y"


#' @examples
#' data(Trentino_hourly_T)
#' date<-date.time(day.begin="01/01/2004", day.end="31/12/2005", date.format= "d/m/y")

#' @seealso \code{\link{Th_int_series}}


####################################################
# GENERATES DATES AND HOURS FROM START TO END
# REQUIRES LIBRARY chron
####################################################

date.time<-function(day.begin, day.end, h.begin=0, h.end=23, date.format= "y/m/d")

{


h.begin.mm_ss<-paste(as.character(h.begin), ":00:00", sep="")
h.end.mm_ss<-paste(as.character(h.end), ":00:00", sep="")
gg.hh_ini<-chron(dates.=day.begin, times.=h.begin.mm_ss, format = c(dates =date.format, times = "h:m:s"))
gg.hh_fin<-chron(dates.=day.end, times.=h.end.mm_ss, format = c(dates = date.format, times = "h:m:s"))

index<-gg.hh_ini+seq(from=0, to=(gg.hh_fin-gg.hh_ini), by=1/24)

year<-as.numeric(as.character(years(index)))
month<-as.numeric(months(index))
day<-as.numeric(days(index))
hour<-hours(index)
date.time<-data.frame(year=year, month=month, day=day, hour=hour)
return(date.time)
}