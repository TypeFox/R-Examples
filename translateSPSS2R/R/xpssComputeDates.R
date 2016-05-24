#' Calculates the difference of time between two dates  in days
#'
#'
#' @description R Implementation of the SPSS \code{CTIME.DAYS} Function. \code{computeCtime_days} is a helper function for xpssCompute. 
#'
#' @usage computeCtime_days(x,date=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @param date atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @details The input format of the date's is "YYYY-MM-DD" as Date, and with timeparameters YYYY-MM-DD HH:MM:SS.
#' 
#' \strong{Examples:} 
#' \tabular{rlll}{
#' \tab Date w/o Time \tab Date with Time 
#' \cr \tab YYYY-MM-DD \tab YYYY-MM-DD HH:MM:SS
#' \cr \tab 2015-01-01 \tab 2015-01-01 00:00:00
#' \cr \tab 2015-12-31 \tab 2015-12-31 23:59:59}
#' 
#' It is possible to calculate the difference between the dates, even if they do not have the same structure. 
#' It is possible to match dates without time paramters, with dates which contain the time parameters. Additionally it is possible to match dates which are not in the right format. 
#' 
#' \strong{Note:} The calculation of two dates, one with the format YYYY-MM-DD and the other one with the format DD-MM-YYYY is valid, but not recommended.
#' 
#' @return \code{computeCtime_days} returns the difference between x and date in days.
#' @author Bastian Wiessner
#' @seealso \code{\link{difftime}} \code{\link{DateTimeClasses}} \code{\link{as.POSIXlt}} \code{\link{strptime}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2013-09-14", fun="computeCtime_days", date="2013-09-11") 
#' xpssCompute(x="2015-11-28", fun="computeCtime_days", date="2014-11-28") 
#' 
#' @export
#' 
#' 

computeCtime_days <- function(x,date=NULL){
  
  #check for right date format
  
  splitdate1 <- str_split(date,pattern = "-")
  if(length(splitdate1[[1]])<3){
    stop("wrong date format at date")
  }
  splitdate2 <- str_split(x,pattern = "-")
  if(length(splitdate2[[1]])<3){
    stop("wrong date format at x")
  }
  if(str_length(splitdate1[[1]][1])==2){
    date <- as.Date(date, "%d-%m-%Y") 
  } else {
    date <- as.Date(strptime(date, "%Y-%m-%d"))
  }
  if(str_length(splitdate2[[1]][1])==2){
    x <- as.Date(x, "%d-%m-%Y") 
  } else {
    x <- as.Date(strptime(x, "%Y-%m-%d"))
  }
  
  out <- difftime(x,date,units="days")

  out <- base::round(out)
  return(out)
}



#' Calculates the difference between two dates in hours
#'
#'
#'  R Implementation of the SPSS \code{CTIME.HOURS} Function. \code{computeCtime_hours} is a helper function for xpssCompute. 
#'
#'
#' @usage computeCtime_hours(x,date)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @param date atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @details The input format of date's with timeparameter's is YYYY-MM-DD HH:MM:SS. All parameters are necessary to calculate the hour difference!
#' 
#' 
#' @return Returns the difference between date and x in hours.
#' @author Bastian Wiessner
#' @seealso \code{\link{difftime}} \code{\link{DateTimeClasses}} \code{\link{as.POSIXlt}} \code{\link{strptime}}
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x="2013-09-14 12:12:12", fun="computeCtime_hours", date="2013-09-14 10:10:10") 
#' xpssCompute(x="2013-09-14 12:12:12", fun="computeCtime_hours", date="2013-09-06 22:10:10") 
#'
#' @export
#'
#'


computeCtime_hours <- function(x,date){
  if(is.na(strptime(date, "%Y-%m-%d %H:%M:%S")) || is.na(strptime(x, "%Y-%m-%d %H:%M:%S"))){
    stop("Date format is not correct. The right format is YYYY-MM-DD HH:MM:SS")
  } else {
    out <-   difftime(strptime(x, format = "%Y-%m-%d %H:%M:%S"),
                      strptime(date, format = "%Y-%m-%d %H:%M:%S"),units="hours")
    out <- base::round(out)  
  }  
  return(out)
}



#' Calculates the difference between two dates in minutes
#'
#'
#' R Implementation of the SPSS \code{CTIME.MINUTES} Function.\code{computeCtime_minutes} is a helper function for xpssCompute.
#'
#'
#' @usage computeCtime_minutes(x,date)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @param date atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @details The input format of date's with timeparameter's is YYYY-MM-DD HH:MM:SS. All parameters are necessary to calculate the minute difference!
#'
#' @return Returns the difference between x and date in minutes.
#' @author Bastian Wiessner
#' @seealso \code{\link{difftime}} \code{\link{DateTimeClasses}} \code{\link{as.POSIXlt}} \code{\link{strptime}}
#' @keywords internal
#' @examples
#'
#' xpssCompute(x="2013-09-14 12:12:12", fun="computeCtime_minutes", date="2013-09-14 10:10:10") 
#' xpssCompute(x="2013-09-14 12:12:12", fun="computeCtime_minutes", date="2013-09-06 22:10:10") 
#'
#' @export


computeCtime_minutes <- function(x,date){
  if(is.na(strptime(date, "%Y-%m-%d %H:%M:%S")) || is.na(strptime(x, "%Y-%m-%d %H:%M:%S"))){
    stop("Date format is not correct. The right format is YYYY-MM-DD HH:MM:SS")
  } else {
    out <-   difftime(strptime(x, format = "%Y-%m-%d %H:%M:%S"),
                      strptime(date, format = "%Y-%m-%d %H:%M:%S"),units="mins")
    out <- base::round(out)    
  }
  return(out)
}

#' Calculates the difference between two dates in seconds
#'
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CTIME.SECONDS} Function. 
#'
#'
#' @usage computeCtime_seconds(x,date)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @param date atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data 
#' @details The input format of date's with timeparameter's is YYYY-MM-DD HH:MM:SS. All parameters are necessary to calculate the second difference!
#'
#' @return Returns the difference between date and x in seconds.
#' @author Bastian Wiessner
#' @seealso \code{\link{difftime}} \code{\link{DateTimeClasses}} \code{\link{as.POSIXlt}} \code{\link{strptime}}
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x="2013-09-14 12:12:12", fun="computeCtime_seconds", date="2013-09-14 10:10:10") 
#' xpssCompute(x="2013-09-14 12:12:12", fun="computeCtime_seconds", date="2013-09-06 22:10:10") 
#'
#' @export


computeCtime_seconds <- function(x,date){
  if(is.na(strptime(date, "%Y-%m-%d %H:%M:%S")) || is.na(strptime(x, "%Y-%m-%d %H:%M:%S"))){
    stop("Date format is not correct. The right format is YYYY-MM-DD HH:MM:SS")
  } else {
    out <-   difftime(strptime(x, format = "%Y-%m-%d %H:%M:%S"),
                      strptime(date, format = "%Y-%m-%d %H:%M:%S"),units="secs")
  out <- base::round(out)
  }
  return(out)
}


#' Creates a  date with the format day-month-year
#'
#'
#' R Implementation of the SPSS \code{DATE.DMY} Function. \code{computeDate_dmy} is a helper function for xpssCompute.
#'
#'
#' @usage computeDate_dmy(day=NULL,month=NULL, year= NULL)
#' @param day atomic numeric or integer.
#' @param month atomic numeric or integer.
#' @param year atomic numeric or integer.

#' @details An character string as Date. Returns a date value corresponding to the indicated day, month, and year. The arguments must resolve to integers or numerics, with day between 1 and 31, month between 1 and 12, and year a four-digit integer value.
#' @return Returns a date object vector of the structure day-month-year.
#'
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_mdy}} \code{\link{computeDate_moyr}}
#' @importFrom stringr str_length 
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x=5, fun="computeDate_dmy", month=12,year=2006) 
#' xpssCompute(x=10, fun="computeDate_dmy", month=10,year=2010) 
#' 
#' @export


computeDate_dmy <- function(day=NULL,month=NULL, year= NULL){
  options(warn = -1)
  if((!(is.numeric(day))) || (!(is.numeric(month))) || (!(is.numeric(year)))){
    stop("day month and year has to be numeric")
  }
  if(length(day) == length(month) && length(month) == length(year)){
    if(day>31 || day<1 || month>12 || month<1)
    {stop("day or month does not exist")}
    if(str_length(round(day))<2){
      day <- paste0(0,round(day))
    }
    if(str_length(month)<2){
      month <- paste0(0,round(month))
    }
    x <- paste0(day,"-",month,"-",year)  
  }
  options(warn = 0)
  return(x)
}


#' Creates a date with the format month-day-year
#'
#'
#' R Implementation of the SPSS \code{DATE.MDY} Function. \code{computeDate_mdy} is a helper function for xpssCompute. 
#'
#'
#' @usage computeDate_mdy(month=NULL,day=NULL, year= NULL)
#' @param day atomic numeric or integer.
#' @param month atomic numeric or integer.
#' @param year atomic numeric or integer.

#' @details An character string as Date. Returns a date value corresponding to the indicated day, month, and year. The arguments must resolve to integers or numerics, with day between 1 and 31, month between 1 and 12, and year a four-digit integer value.
#'
#' @return Returns a character or character vector of the structure month-day-year 
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_dmy}} \code{\link{computeDate_moyr}}
#' @importFrom stringr str_length 
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x=5, fun="computeDate_dmy", day=12,year=2006) 
#' xpssCompute(x=10, fun="computeDate_dmy", day=10,year=2010)

#'
#' @export


computeDate_mdy <- function(month=NULL, day=NULL, year= NULL){
  if((!(is.numeric(day))) || (!(is.numeric(month))) || (!(is.numeric(year)))){
    stop("day month and year has to be numeric")
  }
  if(day>31 || day<1 || month>12 || month<1)
  {stop("day or month does not exist")}
  
  if(str_length(day)<2){
    day <- paste0(0,day)
  }
  if(str_length(month)<2){
    month <- paste0(0,month)
  }
  out <- paste0(month,"-",day,"-",year)
  
  return(out)
}



#' Creates a date with the format month-year
#'
#'
#' R Implementation of the SPSS \code{DATE.MOYR} Function. \code{computeDate_moyr} is a helper function for xpssCompute.  
#'
#'
#' @usage computeDate_moyr(month=NULL, year= NULL)
#'
#' @param month atomic numeric or integer.
#' @param year atomic numeric or integer.

#' @details An character string as Date. Returns a date value corresponding to the indicated month, and year. The arguments must resolve to integers or numerics, with month between 1 and 12, and year a four-digit integer value.
#'
#' @return Returns a character or character vector of the structure month-year. 
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_dmy}} \code{\link{computeDate_mdy}}
#' @importFrom stringr str_length 
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x=5, fun="computeDate_moyr", year=2006) 
#' xpssCompute(x=10, fun="computeDate_moyr",year=2010)
#'
#' @export

computeDate_moyr <- function(month=NULL, year= NULL){
  if((!(is.numeric(month))) || (!(is.numeric(year)))){
    stop("day month and year has to be numeric")
  }
  if(month>12 || month<1)
  {stop("month does not exist")}

  if(str_length(month)<2){
    month <- paste0(0,month)
  }
    out <- paste0(month,"-",year)
  
  return(out)
}


#' Creates a date with the format year/quarter
#'
#'
#' R Implementation of the SPSS \code{DATE.QYR} Function. \code{computeDate_qyr} is a helper function for xpssCompute.
#'
#' @usage computeDate_qyr(day=NULL,month=NULL, year= NULL)
#' @param day atomic numeric or integer.
#' @param month atomic numeric or integer.
#' @param year atomic numeric or integer.

#' @details An character string as Date. Returns a date value corresponding to the indicated month, and year. The arguments must resolve to integers or numerics, with day between 1 and 31, month between 1 and 12, and year a four-digit integer value.
#'
#' @return Returns a atomic character of the structure year / quarter. 
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_wkyr}} \code{\link{computeDate_yrday}}
#' @importFrom zoo as.yearqtr
#' @importFrom stringr str_length 
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x=5, fun="computeDate_qyr", month=10, year=2006) 
#' xpssCompute(x=10, fun="computeDate_qyr", month=1, year=2010)
#'
#' @export



computeDate_qyr <- function(day=NULL,month=NULL, year= NULL){
  if((!(is.numeric(day))) || (!(is.numeric(month))) || (!(is.numeric(year)))){
    stop("day month and year has to be numeric")
  }
  if(day>31 || day<1 || month>12 || month<1)
  {stop("day or month does not exist")}
   
  if(str_length(day)<2){
    day <- paste0(0,day)
  }
  if(str_length(month)<2){
    month <- paste0(0,month)
  }
  
  out <- strptime(paste0(day,"-",month,"-",year),format = "%d-%m-%Y")
  out <- as.yearqtr(out)
  out <- format(x=out,format = "%Y/0%q")
  return(out)
}


#' Creates a date with the format year/calendar week
#'
#'
#' R Implementation of the SPSS \code{DATE.WKYR} Function. \code{computeDate_wkyr} is a helper function for xpssCompute.
#'
#' @usage computeDate_wkyr(day=NULL, month=NULL, year= NULL)
#' @param day atomic numeric or integer.
#' @param month atomic numeric or integer.
#' @param year atomic numeric or integer.

#' @details An character string as Date. Returns a date value corresponding to the indicated month, and year. The arguments must resolve to integers or numerics, with day between 1 and 31, month between 1 and 12, and year a four-digit integer value.
#'
#' @return Returns a atomic character of the structure year / calendar week. 
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_qyr}} \code{\link{computeDate_yrday}}
#' @importFrom stringr str_length
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x=5, fun="computeDate_wkyr", month=10, year=2006) 
#' xpssCompute(x=10, fun="computeDate_wkyr", month=1, year=2010)
#' 
#'
#' @export

computeDate_wkyr <- function(day=NULL, month=NULL, year= NULL){
  if((!(is.numeric(day))) || (!(is.numeric(month))) || (!(is.numeric(year)))){
    stop("week and year has to be numeric")
  }
  if(month>12 || month<1)
  {stop("week does not exist")}
  if(day>32 || day<1)
  {stop("week does not exist")}

  out <- paste0(day,"-",month,"-",year)
  out <- as.Date(out, format="%d-%m-%Y")
  week <- as.numeric( format(out+3, "%U"))
  out <- paste0(format(x=out,format = "%Y/w"),week)
  return(out)
}


#' Creates a date with the format year/yearday
#'
#'
#' R Implementation of the SPSS \code{DATE.YRDAY} Function. \code{computeDate_yrday} is a helper function for xpssCompute.
#'
#' @usage computeDate_yrday(day=NULL,month=NULL, year= NULL)
#' @param day atomic numeric or integer.
#' @param month atomic numeric or integer.
#' @param year atomic numeric or integer.

#' @details An character string as Date. Returns a date value corresponding to the indicated month, and year. The arguments must resolve to integers or numerics, with day between 1 and 31, month between 1 and 12, and year a four-digit integer value.
#'
#' @return Returns a atomic character of the structure year/yearday. 
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_qyr}} \code{\link{computeDate_wkyr}}
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x=5, fun="computeDate_yrday", month=10, year=2006) 
#' xpssCompute(x=10, fun="computeDate_yrday", month=1, year=2010)
#' @export


computeDate_yrday <- function(day=NULL,month=NULL, year= NULL){
  if((!(is.numeric(day))) || (!(is.numeric(month))) || (!(is.numeric(year)))){
    stop("day month and year has to be numeric")
  }
  if(day>31 || day<1 || month>12 || month<1)
  {stop("day or month does not exist")}
  
  if(str_length(day)<2){
    day <- paste0(0,day)
  }
  if(str_length(month)<2){
    month <- paste0(0,month)
  }
  
  out <- paste0(year,"-",month,"-",day)
  doy <- strftime(out, format = "%j")
  out <- paste0(year,"-",doy)
  return(out)
}


#' Calculates the number of passed hours on basis of the given days
#'
#'
#' R Implementation of the SPSS \code{TIME.DAYS} Function. \code{computeTime_days} is a helper function for xpssCompute.
#'
#' @usage computeTime_days(x=NULL)
#' @param x atomic numeric or integer.
#'
#' @return Returns the number of passed hours as atomic numeric. 
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_qyr}} \code{\link{computeDate_wkyr}}
#' @keywords internal
#' @examples
#' xpssCompute(x=5, fun="computeTime_days") 
#' xpssCompute(x=12, fun="computeTime_days")
#'
#' @export

computeTime_days <- function(x=NULL){
  if(is.numeric(x)){
    out <- x*24  
  } else {
    stop("x argument has to be numeric")
  }
  
  return(out)
}


#' Creates a date with the format Hour-Minute-Second
#'
#'
#' R Implementation of the SPSS \code{TIME.HMS} Function. \code{computeTime_hms} is a helper function for xpssCompute.
#'
#' @usage computeTime_hms(hour=NULL,minute=NULL,second=NULL)
#' @param hour atomic numeric or integer.
#' @param minute atomic numeric or integer.
#' @param second atomic numeric or integer.
#'
#' @return Returns a character string as Timedate.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeDate_qyr}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_length
#' @keywords internal
#' @examples
#' 
#' xpssCompute(x=4, fun="computeTime_hms", minute= 10, second=23) 
#' xpssCompute(x=14, fun="computeTime_hms", minute= 34, second=1) 
#'
#' @export

computeTime_hms <- function(hour=NULL,minute=NULL,second=NULL){
  if((!(is.numeric(hour))) || (!(is.numeric(minute))) || (!(is.numeric(second)))){
    stop("hour, minute and second has to be numeric") 
  }
   
  if(hour>24 || hour<0){
    stop("hour value is limit between 0 and 24")
  }
  if(minute>60 || minute<0){
    stop("minute value is limit between 0 and 24")
  }
  if(second>60 || second<0){
    stop("second values is limit between 0 and 24")
  }
  
  if(str_length(round(hour))<2){
    hour <- paste0(0,round(hour))
  }
  if(str_length(minute)<2){
    minute <- paste0(0,round(minute))
  }
  if(str_length(second)<2){
    second <- paste0(0,round(second))
  }
  
  out <- paste0(hour,":",minute,":",second)
  return(out)
}


#' Extracts the date out of a date string.
#'
#'
#' R Implementation of the SPSS \code{XDATE.DATE} Function. \code{computeXdate_date} is a helper function for xpssCompute.
#'
#' @usage computeXdate_date(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_date extract the complete date out of the date string, the time componentes will be omitted.
#' @return Returns a character string with date values.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_hour}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_length str_split
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:22:20", fun="computeXdate_date") 
#' xpssCompute(x="2022-02-20 21:22:12", fun="computeXdate_date")
#' 
#' @export

computeXdate_date <- function(x = NULL){ 
  
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- as.Date(strptime(x, "%d-%m-%Y"))
  } else {
    x <- as.Date(x, "%Y-%m-%d") 
  }
  if(is.null(x)){
    stop("Maybe the input has the wrong format?")
  }
  return(x)
}

#' Extracts the hour value out of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.HOUR} Function. \code{computeXdate_hour} is a helper function for xpssCompute.
#'
#' @usage computeXdate_hour(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_hour extract the hour component out of the date string, the other date and time componentes will be omitted.
#' @return Returns a character string with the hour value values.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:22:20", fun="computeXdate_hour") 
#' xpssCompute(x="2022-02-20 21:22:12", fun="computeXdate_hour")
#'
#' @export

computeXdate_hour <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format at date")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y %H:%M:%S")
  } else {
    x <- strptime(x, "%Y-%m-%d %H:%M:%S") 
  }
  out <- format(x, "%H") 
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}



#' Calculates the date of year on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.JDAY} Function. \code{computeXdate_jday} is a helper function for xpssCompute.
#'
#' @usage computeXdate_jday(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_jday caluclates the date of year out of the date string, the other date and time componentes will be omitted.
#' @return Returns a character string with the date of year.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:22:20", fun="computeXdate_jday") 
#' xpssCompute(x="2022-02-20 21:22:12", fun="computeXdate_jday")
#'
#'
#' @export

computeXdate_jday <- function(x = NULL){
  #' @importFrom stringr str_split str_length
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y")
  } else {
    x <- strptime(x, "%Y-%m-%d") 
  }  
  out <- format(x, "%j")
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}



#' Extracts the date of month on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.MDAY} Function. \code{computeXdate_mday} is a helper function for xpssCompute.
#'
#' @usage computeXdate_mday(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_mday extract the day component out of the date string, the other date and time componentes will be omitted.
#' @return Returns a character string with the date's day.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:22:20", fun="computeXdate_mday") 
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_mday")
#'
#' @export

computeXdate_mday <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y")
  } else {
    x <- strptime(x, "%Y-%m-%d") 
  } 
  
  out <- format(x, "%d")
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}



#' Extracts the minute component on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.MINUTE} Function. \code{computeXdate_minute} is a helper function for xpssCompute.
#'
#' @usage computeXdate_minute(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_minute extract the minute component out of the date string, the other date and time componentes will be omitted.
#' @return Returns a character string with the minute component.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_minute") 
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_minute")

#' @export

computeXdate_minute <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }  
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y %H:%M:%S")
  } else {
    x <- strptime(x, "%Y-%m-%d %H:%M:%S") 
  }
    out <- as.numeric(format(x, "%M"))
  if(is.null(out)){
    stop("Minute component is missing. Maybe the input has the wrong format?")
  }
  return(out)
}


#' Extracts the month component on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.MINUTE} Function. \code{computeXdate_month} is a helper function for xpssCompute.
#'
#' @usage computeXdate_month(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_month extract the month component out of the date string, the other date and time componentes will be omitted.
#' @return Returns a character string with the month component.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_month") 
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_month")
#'
#' @export

computeXdate_month <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y")
  } else {
    x <- strptime(x, "%Y-%m-%d") 
  } 
  out <- as.numeric(format(x, "%m"))
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}


#' Calculates the quarter on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.QUARTER} Function. \code{computeXdate_quarter} is a helper function for xpssCompute.
#'
#' @usage computeXdate_quarter(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_quarter calculates the quarter component out of the date string, the other date and time componentes will be omitted.
#' @return Returns a character string with the quarter component.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length str_extract
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_quarter") 
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_quarter")
#' @export

computeXdate_quarter <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y")
  } else {
    x <- strptime(x, "%Y-%m-%d") 
  } 
  out <- quarters(x=x)
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  out <- as.numeric(str_extract(out, "[1-4]"))
  return(out)
}



#' Extracts the second on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.SECOND} Function. \code{computeXdate_second} is a helper function for xpssCompute.
#'
#' @usage computeXdate_second(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_second calculates the second component out of the date string, the other date and time componentes will be omitted.
#' @return Returns a character string with the second component.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_second") 
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_second")
#'
#' @export

computeXdate_second <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y %H:%M:%S")
  } else {
    x <- strptime(x, "%Y-%m-%d %H:%M:%S") 
  }
  out <- format(x, "%S")
  if(is.null(out)){
    stop("Second component is missing. Maybe the input has the wrong format?")
  }
  return(out)
}



#' Calculates the difference of days between the entered date and October 14, 1582.
#'
#'
#' R Implementation of the SPSS \code{XDATE.DAY} Function. \code{computeXdate_tday} is a helper function for xpssCompute.
#'
#' @usage computeXdate_tday(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_tday calculates the difference of days between October 14, 1582 and the entered date.
#' @return Returns a character string with difference of days between the dates.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_tday") 
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_tday")

#'
#' @export

computeXdate_tday <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y %H:%M:%S")
  } else {
    x <- strptime(x, "%Y-%m-%d %H:%M:%S") 
  } 
  
  gregorian <- as.POSIXct(x = "1582-10-14",tz = "")
    
  out <- difftime(time1 = x,time2 = gregorian,units = "days")
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}

#' Extracts the time componente on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.TIME} Function. \code{computeXdate_time} is a helper function for xpssCompute.
#'
#' @usage computeXdate_time(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_time extract the time components out of the date string, the other date components will be omitted.
#' @return Returns a character string with the time components.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_time") 
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_time")
#' 
#' @export

computeXdate_time <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y %H:%M:%S")
  } else {
    x <- strptime(x, "%Y-%m-%d %H:%M:%S") 
  } 
  
  out <- format(x, "%H:%M:%S")
  if(is.null(out)){
    stop("Time components are missing. Maybe the input has the wrong format?")
  }
  return(out)
}

#' Calcualtes the calendar week  on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.WEEK} Function. \code{computeXdate_week} is a helper function for xpssCompute.
#'
#' @usage computeXdate_week(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_week caluclates the calendar week on basis of the given date string.
#' @return Returns a character string with the calendar week.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeDate_wkyr}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_week") 
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_week")
#'
#' @export

computeXdate_week <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format at date")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y")
  } else {
    x <- strptime(x, "%Y-%m-%d") 
  } 
  
  out <- format(x, "%W")
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}

#' Calcualtes the day of week on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.WEEK} Function. \code{computeXdate_wkday} is a helper function for xpssCompute.
#'
#' @usage  computeXdate_wkday(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_wkday caluclates the calendar week on basis of the given date string. Result is a number between 0 and 6. 0 stands for Sunday, 6 for Saturday.
#' @return Returns a character string with the day of week.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeXdate_wkday}} \code{\link{computeXdate_year}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20", fun="computeXdate_wkday")
#' xpssCompute(x="2022-10-20 21:22:12", fun="computeXdate_wkday")
#'
#' @export

computeXdate_wkday <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y")
  } else {
    x <- strptime(x, "%Y-%m-%d") 
  } 
  
  out <- format(x, "%w")
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}

#' Extracts the year on basis of a given date
#'
#'
#' R Implementation of the SPSS \code{XDATE.YEAR} Function. \code{computeXdate_year} is a helper function for xpssCompute.
#'
#' @usage computeXdate_year(x=NULL)
#' @param x atomic object of class \code{character, POSIXlt  or  POSIXt} holding date data.
#' @details computeXdate_year extracts the year on basis of the given date string. 
#' @return Returns a character string with the day of week.
#' @author Bastian Wiessner
#' @seealso \code{\link{computeXdate_date}} \code{\link{computeXdate_wkday}}
#' @importFrom stringr str_split str_length
#' @keywords internal
#' @examples
#' xpssCompute(x="2015-02-15 20:11:20",fun="computeXdate_year")
#' xpssCompute(x="2022-10-20 21:22:12",fun="computeXdate_year")
#'
#' @export


computeXdate_year <- function(x = NULL){
  splitdate <- str_split(x,pattern = "-")
  if(length(splitdate[[1]])<3){
    stop("wrong date format")
  }
  if(str_length(splitdate[[1]][1])==2){
    x <- strptime(x, "%d-%m-%Y")
  } else {
    x <- strptime(x, "%Y-%m-%d") 
  } 
  
  out <- format(x, "%Y")
  if(is.null(out)){
    stop("Maybe the input has the wrong format?")
  }
  return(out)
}
