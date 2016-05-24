#######################################################################################################
#               Written by Andrew R Sommerlot <andrewrs@vt.edu>, March 2015                           #
#######################################################################################################
#' Converts a date time series to SWAT IO format dates
#' @param dates - Input vector of dates in as Date or character class.
#' @param dateformat - Format of input dates. Can be 'yyyy-mm-dd', 'mm-dd-yyyy', 'mm/dd/yyyy', 'yyyy/mm/dd','dd/mm/yyy', or 'dd-mm-yyy'. Default is 'yyyy-mm-dd'.
#' @return returns a time series of swat IO format dates from the input
#' @examples
#' dates = c('2000-12-28', '2000-12-29')
#' getSWATdates(dates)
#' @export


getSWATdates <- function(dates, dateformat = 'yyyy-mm-dd'){

  if(dateformat == 'yyyy-mm-dd'){
    dateformat <- c('year', 'month', 'day')
    colsep <- '-'
  } else if(format == 'mm-dd-yyyy'){
    dateformat <- c('month', 'day', 'year')
    colsep <- '-'
  } else if(format == 'mm/dd/yyyy'){
    dateformat <- c('month', 'day', 'year')
    colsep <- '/'
  } else if(format == 'yyyy/mm/dd'){
    dateformat <- c('year', 'month', 'day')
    colsep <- '/'
  } else if(format == 'dd/mm/yyy'){
    dateformat <- c('day', 'month', 'year')
    colsep <- '/'
  }else if(format == 'dd-mm-yyy'){
    dateformat <- c('day', 'month', 'year')
    colsep <- '/'
  }else {
    stop(c("Bad Date format", "use 'yyyy-mm-dd', 'yyyy/mm/dd', 'mm-dd-yyyy', 'dd/mm/yyy', 'dd-mm-yyy' or 'mm-dd-yyyy\n'", "Make sure input matches one of these formats and 'dateformat' is set accordingly"))
  }


  ####added from original code to satisfy the colsplit requirement of character object. #remember as.___ does NOT coerce, it needs a to be set to a object
  dates1 <- as.character(dates)


  ###############################################################################################
  ###############################################################################################

  datesm <- data.frame(matrix(nrow = length(dates1), ncol = 3))

  datessplit <- strsplit(dates1, colsep)

  for(i in 1:length(dates1)){
    datesm[i,1] <- datessplit[[i]][1]
    datesm[i,2] <- datessplit[[i]][2]
    datesm[i,3] <- datessplit[[i]][3]

  }

  dates1 <- data.frame(datesm, stringsAsFactors = FALSE)

  dates1[,1] <- as.numeric(dates1[,1])
  dates1[,2] <- as.numeric(dates1[,2])
  dates1[,3] <- as.numeric(dates1[,3])

  dates1 <- as.data.frame(dates1)

  colnames(dates1) <- c(dateformat)

  ###############################################################################################
  ###############################################################################################
  #Create the 2nd to last data frame for dates
  dates.format <- data.frame(matrix(ncol = 2, nrow = nrow(dates1)))
  colnames(dates.format) <- c('year', 'days')

  #create the add days matrix. this defines how many extra days are added given the month
  add.days <- data.frame(matrix(ncol = 2, nrow = 12))
  colnames(add.days) <- c('month', 'days')

  #define the months and corresponing extra days to add. ex.days is calculated by placing the amount of days from the ith-1 month.
  ex.days <- c(0,31,59,90,120,151,181,212,243,273,304,334)
  months <- c(1,2,3,4,5,6,7,8,9,10,11,12)

  #To take care of the leap year change:
  leap.years <- vector()
  leap.years <- as.numeric(leap.years)
  year <- 1900
  for(i in 1:100){
    year <- year+ (1*4)
    leap.years[i] <- year
  }

  ln.leap.years <- length(leap.years)

  for(i in 1:ln.leap.years){
    if(dates1$year[1] == leap.years[i]){
      ex.days[3] <- 29
    }
  }

  #fill in add.days
  add.days$month <- months
  add.days$days <- ex.days



  extra.days <- vector()
  extra.days <- as.numeric(extra.days)
  ex.day <- 1

  #define extra days and fill in the days in dates.format
  for(i in 1:nrow(dates1)){
    month <- dates1$month[i]
    for(j in 1:12){
      if(month == add.days$month[j]){
        ex.day <- add.days$days[j]
      }
      extra.days[i] <- ex.day
    }
  }


  #fill out the day column of dates.format
  dates.format$days <- dates1$day + extra.days
  dates.format$year <- dates1$year

  #format the day column to '000' format
  dates.format$days <- formatC(dates.format$days, width = 3, flag = '0')

  #create new column in dates.format that is equivalent to the SWAT date
  dates.format$fullDate <- paste(dates.format$year, dates.format$days, sep = "")
  date.full <- dates.format$fullDate

  return(date.full)

}




