#######################################################################################################
#               Written by Andrew R Sommerlot <andrewrs@vt.edu>, Februrary 2015                       #
#######################################################################################################
#' generates a vector of dates in the swat model IO format
#' @param startDate - Start date for output vector, supports common date formats
#' @param endDate - End date for output vector, supports common date formats
#' @return returns a generated time series of swat IO format dates
#' @examples
#' startDate = '2000-12-28'
#' endDate = '2000-12-29'
#' genSWATdates(startDate=startDate, endDate=endDate)
#' @export



genSWATdates <- function(startDate, endDate){


  #get the single date
  date <- as.Date(startDate)

  #calculate the length of the generated date file
  edate <- as.Date(endDate)

  dates <- seq.Date(date, edate, by <- 'day')

  dates <- as.character(dates)

  datesm <- data.frame(matrix(nrow = length(dates), ncol <- 3))

  datessplit <- strsplit(dates, '-')

  for(i in 1:length(dates)){
    datesm[i,1] <- as.numeric(datessplit[[i]][1])
    datesm[i,2] <- as.numeric(datessplit[[i]][2])
    datesm[i,3] <- as.numeric(datessplit[[i]][3])

  }

  dates <- datesm

  colnames(dates) <- c('year', 'month', 'day')


  #Create the 2nd to last data frame for dates
  dates.format <- data.frame(matrix(ncol = 2, nrow = nrow(dates)))
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
  for(j in 1:900){
    year <- year+ (1*4)
    leap.years[j] <- year
  }

  ln.leap.years <- length(leap.years)

  for(k in 1:ln.leap.years){
    if(dates$year[1] == leap.years[k]){
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
  for(l in 1:nrow(dates)){
    month <- dates$month[l]
    for(m in 1:12){
      if(month == add.days$month[m]){
        ex.day <- add.days$days[m]
      }
      extra.days[l] <- ex.day
    }
  }


  #fill out the day column of dates.format
  dates.format$days <- dates$day + extra.days
  dates.format$year <- dates$year

  #format the day column to '000' format
  dates.format$days <- formatC(dates.format$days, width = 3, flag = '0')

  #create new column in dates.format that is equivalent to the SWAT date
  dates.format$fullDate <- paste(dates.format$year, dates.format$days, sep = "")
  date.full <- dates.format$fullDate

  return(date.full)

}
