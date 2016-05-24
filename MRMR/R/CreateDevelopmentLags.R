CreateLagsFromIntegers = function(LagValues, DevelopmentPeriod = months(1), Verbose = TRUE)
{
  if (!is.period(DevelopmentPeriod)) stop ("A period object was not specified for the DevelopmentPeriod.")
  
  if (!is.integer(LagValues)){
    if (Verbose) warning ("LagValues converted to integer.")
    LagValues = as.integer(LagValues)
  }
  
  DevelopmentLag = LagValues * DevelopmentPeriod
  
}

CreateLagsFromEvalDates = function(EvalDates, OriginPeriods, DevelopmentPeriod = months(1))
{
  if (!is.period(DevelopmentPeriod)) stop ("A period object was not specified for the DevelopmentPeriod.")
  
  if (is.null(OriginPeriods)) stop ("OriginPeriods were not specified.")
  
  DevelopmentInterval = new_interval(int_start(OriginPeriods), EvalDates + days(3))
  DevelopmentLag = suppressMessages (DevelopmentInterval / DevelopmentPeriod)
  DevelopmengLag = DevelopmentLag * DevelopmentPeriod
}

#' Create triangle development lags
#' 
#' If the triangle dataframe does not record development lags as lubridate periods, they must be 
#' created. 
#' Development lags may be established one of three ways:
#' 1. The development lags are passed in as lubridate periods.
#'    Everything's cool. The evaluation dates are established by adding the periods to the starting 
#'    point
#'    of the origin periods.
#' 2. The development lags are passed in as integers, with a presumed time period.
#'    The program will establish lubridate period objects using the integers and time periods and then 
#'    proceed 
#'    as above.
#' 3. An evaluation date is passed in. Here we must take the difference between the evaluation dates
#'    and the origin periods. We will assume months as the default period. The user may pass in another.
#' @export CreateDevelopmentLags
#' 
#' @param LagValues Vector of development lags
#' @param DevelopmentPeriod A period object indicating the common time period between evaluations
#' @param EvaluationDates A vector of evaluation dates
#' @param OriginPeriods A vector of origin periods
#' @param Verbose Show warnings?
#' 
#' @return A vector of intervals
#' 
#' @seealso \code{\link{CreateDevelopmentLags}}, \code{\link{CreateEvaluationDates}}
#' @examples
#' library(lubridate)
#' 
#' Sys.setenv(TZ='UTC')
#' # Case 1
#' DevelopmentLag = c(months(12), months(24), months(12))
#' DevelopmentLag
#' 
#' # Case 2
#' LagValues = c(12, 24, 12)
#' dPeriod = months(1)
#' DevelopmentLags = CreateDevelopmentLags(LagValues, DevelopmentPeriod = dPeriod)
#' DevelopmentLags
#' 
#' # Case 3
#' OriginStart = c(mdy("1/1/2000"), mdy("1/1/2000"), mdy("1/1/2001"))
#' OriginEnd = c(mdy("12/31/2000"), mdy("12/31/2000"), mdy("12/31/2001"))
#' OriginPeriods = CreateOriginPeriods(OriginStart, OriginEnd)
#' 
#' EvaluationDates = c(mdy("12/31/2000"), mdy("12/31/2001"), mdy("12/31/2001"))
#' DevelopmentLags = CreateDevelopmentLags(DevelopmentPeriod = months(1)
#'                                        , EvaluationDates = EvaluationDates
#'                                        , OriginPeriods = OriginPeriods)
#' DevelopmentLags
#' 
#' DevelopmentPeriod = years(1)
#' DevelopmentLags = CreateDevelopmentLags(DevelopmentPeriod = months(1)
#'                                        , EvaluationDates = EvaluationDates
#'                                        , OriginPeriods = OriginPeriods)
#' DevelopmentLag
CreateDevelopmentLags = function(LagValues, DevelopmentPeriod = months(1), EvaluationDates = NULL, OriginPeriods = NULL, Verbose = TRUE)
{
  if (!is.period(DevelopmentPeriod)) stop ("A period object was not specified for the DevelopmentPeriod.")
  
  if (is.null(EvaluationDates))
  {
    return (CreateLagsFromIntegers(LagValues, DevelopmentPeriod, Verbose))
  } else {
    return (CreateLagsFromEvalDates(EvaluationDates, OriginPeriods, DevelopmentPeriod))
  }
  
}