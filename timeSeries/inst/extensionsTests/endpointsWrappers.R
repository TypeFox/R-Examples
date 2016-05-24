

require(timeSeries)


###############################################################################
# FUNCTION:
#  timeNdayInWeek
#  timeLastBizdayInWeek
# FUNCTION:
#  timeLastDayInMonth
#  timeLastNdayInMonth
#  timeLastBizdayInMonth
#  timeNthNdayInMonth
# FUNCTION:
#  timeLastDayInQuarter
#  timeLastNdayInQuarter
#  timeLastBizdayInQuarter
#  timeNthNdayInQuarter
###############################################################################


###############################################################################
# endpoints
#   extract index values of a given time Series object corresponding to  
#   the last calendarical observation in the specified period 


require(timeSeries)

# Daily and Monthly Series in 2011:
tD <- timeCurrentYear(2011)
tM <- timeCalendar(2011)


###############################################################################
# Weekly Endpoints


# -----------------------------------------------------------------------------
# On Given nDay of Week:


timeNdayInWeek <- function(x, nday=5) {
  X <- align(x)
  DOW <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
  X[dayOfWeek(X) == DOW[nday+1], ]
}


LastDayInWeek <- timeNdayInWeek(tD)
LastDayInWeek
dayOfWeek(LastDayInWeek)


# -----------------------------------------------------------------------------
# Last Bizday In Week:


timeLastBizdayInWeek <- 
  function(x, holidays=holidayNYSE)
{
  # Extend time Sequence:
  x <- timeSequence(
    from = timeDate(x[1]), to = timeDate(x[length(x)]), by = "day")
  
  # Bizdays Function:
  FUN <- function(x, holidays = holidays) {
    holidays <- holidayNYSE()
    posix <- as.POSIXct(x, zone = "", origin = "1970-01-01")
    check <- isBizday(as.timeDate(posix), holidays = holidays, wday = 1:5)
    ans <- rev(x[check])[1]
    ans }
  
  # Create Periods from:
  by <- timeDayInWeek(x, nday=5)
  bySec <- as.numeric(by, "sec")
  xSec <- as.numeric(x, "sec")
  
  # Compute Index:
  INDEX <- findInterval(xSec, bySec + 1)
  INDEX <- INDEX + 1
  is.na(INDEX) <- !(INDEX <= length(by))
  dates <- matrix(apply(matrix(xSec, ncol=1), 2, tapply, INDEX, FUN), ncol=1)
  dates <- as.timeDate(as.POSIXct(dates, zone="GMT", origin="1970-01-01"))
  
  # Return Value:
  dates
}


LastNYBizdayInWeek <- timeLastBizdayInWeek(tD, holidays=holidayNYSE)
LastNYBizdayInWeek
dayOfWeek(LastNYBizdayInWeek)


# Can we attribute the day of week?
ans <- timeSeries(
  data = rnorm(length(LastBizdayInWeek)),
  charvec = timeLastBizdayInWeek(tD),
  recordIDs = data.frame(DOW=dayOfWeek(LastBizdayInWeek)))


###############################################################################
# Monthly Endpoints 


# -----------------------------------------------------------------------------
# Last Calendar Day in Month:
LastDayInMonth <- timeLastDayInMonth(tX, unique=TRUE)
LastDayInMonth
dayOfWeek(LastDayInMonth)


# -----------------------------------------------------------------------------
# Last Friday in Month:
LastFridayInMonth <- timeLastNdayInMonth(tX, nday=5, unique=TRUE)
LastFridayInMonth
dayOfWeek(LastFridayInMonth)


# -----------------------------------------------------------------------------
# Last New-York Bizday in Month
holidayNYSE(2011)
LastNYBizdayInMonth <- timeLastBizdayInMonth(tD, holidays=holidayNYSE(), 
  unique=TRUE)
LastNYBizdayInMonth
dayOfWeek(LastNYBizdayInMonth)


# -----------------------------------------------------------------------------
# 2nd Tuesday in Month:
SecondTuesdayInMonth <- timeNthNdayInMonth(tX, nth=2, nday=2, unique=TRUE)
SecondTuesdayInMonth
dayOfWeek(SecondTuesdayInMonth)


###############################################################################
# Quarterly Endpoints


# -----------------------------------------------------------------------------
# Last Day in Quarter:
timeLastDayInQuarter <- 
  function(charvec, format="%Y-%m-%d", zone="", FinCenter="", unique = FALSE)
{
  # Description:
  #   Returns Last Nday in Quarter
    
  ans <- timeLastDayInMonth(charvec, format, zone = "", FinCenter, unique)
  INDEX <- which ( atoms(ans)[, 2] %in% c(3, 6, 9, 12) )
  ans[INDEX]
}
LastDayInQuarter <- timeLastDayInQuarter(tD, unique=TRUE)  
LastDayInQuarter
dayOfWeek(LastDayInQuarter)


# -----------------------------------------------------------------------------
# Last Friday in Quarter:
timeLastNdayInQuarter <- 
  function(charvec, nday = 1, format = "%Y-%m-%d", zone = "", FinCenter = "", 
    unique = FALSE) 
{
  # Description:
  #   Returns Last Nday first/mid/last inMonths of Quarters
  
  ans <- timeLastNdayInMonth(charvec, nday, format, zone, FinCenter,unique) 
  INDEX <- which ( atoms(ans)[, 2] %in% c(3, 6, 9, 12) )
  ans[INDEX]
}

LastFridayInQuarter <- timeLastNdayInQuarter(tD, nday=5, unique=TRUE) 
LastFridayInQuarter
dayOfWeek(LastFridayInQuarter)


# -----------------------------------------------------------------------------
# Last Bizday in Quarter:
timeLastBizdayInQuarter <- 
  function(charvec, holidays = holidayNYSE(), format = "%Y-%m-%d", 
    zone = "", FinCenter = "", unique = FALSE)
{
  ans <- timeLastBizdayInMonth(charvec, holidays, format, zone, FinCenter, 
    unique)
  INDEX <- which ( atoms(ans)[, 2] %in% c(3, 6, 9, 12) )
  ans[INDEX]
}

LastNYBizdayInQuarter <- timeLastBizdayInQuarter(tD, holidayNYSE(), unique=TRUE)
LastNYBizdayInQuarter
dayOfWeek(LastNYBizdayInQuarter)



# -----------------------------------------------------------------------------
# nth-of Mar/Jun/Sep/Dec Nday in Quarter:
timeNthNdayInQuarter <- 
  function(charvec, nday = 1, nth = 1, inMonths = c(3, 6, 9, 12),
    format = "%Y-%m-%d", zone = "", FinCenter = "", unique = FALSE)
{
  ans <- timeNthNdayInMonth(charvec, nday, nth, format, zone, FinCenter, 
    unique)
  INDEX <- which ( atoms(ans)[, 2] %in% inMonths)
  ans[INDEX]
}


# 2nd Tuesday in (last) Months 3/6/9/12:
SecondTuesdayInLastQuarterMonth <- timeNthNdayInQuarter(tD, nth=2, nday=2, 
  unique=TRUE)
SecondTuesdayInLastQuarterMonth
dayOfWeek(SecondTuesdayInLastQuarterMonth)


# 2nd Tuesday in (first) Months 1/4/7/10:
SecondTuesdayInFirstQuarterMonth <- timeNthNdayInQuarter(tD, nth=2, nday=2, 
  inMonths=c(1, 3, 7, 10), unique=TRUE)
SecondTuesdayInFirstQuarterMonth
dayOfWeek(SecondTuesdayInFirstQuarterMonth)


# IMM Dates:
# The dates are the third Wednesday of March, June, September and December
datesIMM <- timeNthNdayInQuarter(tD, nth=3, nday=3, 
  inMonths=c(3, 6, 9, 12), unique=TRUE)
datesIMM
dayOfWeek(datesIMM)


###############################################################################


