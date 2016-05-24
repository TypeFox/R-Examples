#' Create triangle evaluation dates
#' 
#' Triangle evaluation dates are established by adding development lags to the starting point of the origin periods.
#' @export CreateEvaluationDates
#' 
#' @param OriginPeriod A vector of interval objects
#' @param DevelopmentLag A vector of period objects
#' 
#' @return A vector of intervals
#' 
#' @seealso \code{\link{CreateDevelopmentLags}}, \code{\link{CreateOriginPeriods}}
#' 
#' @examples
#' library(lubridate)
#'
#' OriginStart = c(mdy("1/1/2000"), mdy("1/1/2000"), mdy("1/1/2001"))
#' OriginEnd = c(mdy("12/31/2000"), mdy("12/31/2000"), mdy("12/31/2001"))
#' OriginPeriod = CreateOriginPeriods(OriginStart, OriginEnd) 
#' DevelopmentLag = c(months(12), months(24), months(12))
#' 
#' EvaluationDates = CreateEvaluationDates(OriginPeriod, DevelopmentLag)
#' EvaluationDates

CreateEvaluationDates = function(OriginPeriod, DevelopmentLag)
{
  if (!is.interval(OriginPeriod)) stop ("OriginPeriod is not a valid interval object.")
  if (!is.period(DevelopmentLag)) stop ("DevelopmentLag is not a valid Period object.")
  
  EvaluationDate = int_start(OriginPeriod) + DevelopmentLag - days(1)
}