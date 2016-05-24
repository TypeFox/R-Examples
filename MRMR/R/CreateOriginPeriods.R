OriginPeriodFromDates = function(OriginStart, OriginEnd, Verbose = FALSE)
{
  if(is.null(OriginStart)) stop ("No OriginStart argument supplied.")
  
  if (!is.POSIXt(OriginStart))
    stop ("OriginStart and OriginEnd were supplied, but OriginStart is not a POSIXct or POSIXt class variable.")
  
  if (!is.POSIXt(OriginEnd))
    stop ("OriginStart and OriginEnd were supplied, but OriginEnd is not a POSIXct or POSIXt class variable.")
  
  if (length(OriginStart) != length(OriginEnd)) stop("OriginStart and OriginEnd have different lengths.")
  
  if (sum(OriginStart >= OriginEnd)) stop("At least one OriginStart date comes after OriginEnd.")
  
  OriginPeriod = new_interval(OriginStart, OriginEnd) 
}

OriginPeriodFromYears = function(OriginStartYears, OriginLength = years(1), StartDay = 1, StartMonth = 1, Verbose = FALSE)
{
  if (is.null(OriginStartYears)) stop ("No OriginStartYears argument supplied.")
  
  #  if (is.null(OriginLength)) OriginLength = years(1)
  
  if (class(OriginLength) != "Period") stop ("A Period object was not supplied with the OriginLength argument.")
  
  if (Verbose) warning("You've not supplied a date for the OriginStart. I'm going to assume that you've given me years.")
  
  if (is.double(OriginStartYears)) 
  {
    if (Verbose) warning ("Converting OriginStartYears from double to integer")
    OriginStartYears = as.integer(OriginStartYears)
  }
  
  OriginStart = rep(now(), length(OriginStartYears))
  year(OriginStart) = OriginStartYears
  month(OriginStart) = StartMonth
  day(OriginStart) = StartDay
  hour(OriginStart) = 0
  minute(OriginStart) = 0
  second(OriginStart) = 0
  OriginEnd = OriginStart + OriginLength - days(1)
  OriginPeriod = new_interval(OriginStart, OriginEnd)
  
}

#' @title CreateOriginPeriods
#' 
#' @description
#' This will create a set of origin period values
#' 
#' @details
#' If the triangle dataframe does not have lubridate intervals, they must be created. Origin 
#' periods may be established one of three ways:
#' 1. The origin periods are passed in as POSIX dates.
#'    This is a simple matter of forming the interval using lubridate.
#' 2. The origin periods are passed in with a start date, but no end date.
#'    We need to have a common period to establish the end date.
#' 3. The origin periods are passed in as parts of a date.
#'    This will typically happen if we know the year, but not the month or day. In this case
#'    , the user may pass 
#'    in month and day values
#' @export CreateOriginPeriods
#' 
#' @param OriginStart Either a vector of date-time objects, or a vector of numbers indicating the 
#' year.
#' @param OriginEnd A vector of date-time objects. If this argument is supplied, it is assumed that 
#' OriginStart contains date-time objects.
#' @param OriginLength A Period object. These are easily created as shown in the example below. The 
#' default is a period of one year. If OriginStart and OriginEnd are supplied, this argument is ignored.
#' @param StartDay If OriginStart and OriginEnd are supplied, this argument is ignored.
#' @param StartMonth If OriginStart and OriginEnd are supplied, this argument is ignored.
#' @param Verbose Boolean indicating whether or not to display warning messages.
#' 
#' @return A vector of intervals
#' 
#' @seealso \code{\link{CreateDevelopmentLags}}, \code{\link{CreateEvaluationDates}}
#' 
#' @examples
#'
#' # Case 1
#' library(lubridate)
#' OriginStart = c(mdy("1/1/2000"), mdy("1/1/2000"), mdy("1/1/2001"))
#' OriginEnd = c(mdy("12/31/2000"), mdy("12/31/2000"), mdy("12/31/2001"))
#' 
#' OriginPeriods = CreateOriginPeriods(OriginStart, OriginEnd)
#' OriginPeriods
#' 
#' # Case 2
#' OriginStart = c(mdy("1/1/2000"), mdy("1/1/2000"), mdy("1/1/2001"))
#' OriginPeriods = CreateOriginPeriods(OriginStart, OriginLength = months(12))
#' OriginPeriods
#' 
#' # Case 3
#' OriginStartYear = c(2000, 2000, 2001)
#' OriginPeriods = CreateOriginPeriods(OriginStartYear, OriginLength = years(1)
#'                                      , StartDay = 1, StartMonth = 1)
#' OriginPeriods
#' 
CreateOriginPeriods = function(OriginStart, OriginEnd = NULL
                               , OriginLength = years(1), StartDay = 1
                               , StartMonth = 1, Verbose = FALSE)
{
  if(is.null(OriginStart))
  {
    stop ("No OriginStart argument supplied.")
  }
  
  # Case 1: We've been given dates for the start and end
  if (!is.null(OriginEnd))
  {
    return (OriginPeriodFromDates(OriginStart, OriginEnd))
  }
  
  # Case 2: We've been given a start date, but no end date
  if (is.POSIXt(OriginStart))
  {
    OriginEnd = OriginStart + OriginLength - days(1)
    return (OriginPeriodFromDates(OriginStart, OriginEnd))
  }
  
  # Case 3: We've been given a set of numbers, but no date
  OriginPeriodFromYears(OriginStart, OriginLength, StartDay = 1, StartMonth = 1, Verbose)
}