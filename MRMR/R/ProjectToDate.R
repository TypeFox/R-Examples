#' ProjectToDate
#' 
#' @include TriangleModel.R
#' 
#' @export
#' 
#' @description
#' This function 
#' 
#' @param objTriangleModel A TriangleModel object
#' @param lOriginYears A list of origin years
#' @param AsOfDate A date to which to project
#' 
#' @return
#' A data frame which has projected dates and columns for the new stochastic values
ProjectToDate = function(objTriangleModel, lOriginYears, AsOfDate)
{
  objTriangle = objTriangleModel@Triangle
  DevelopmentInterval = objTriangle@DevelopmentInterval
  
  mojo = lapply(lOriginYears, function(x){
    tz(AsOfDate) = tz(x$EvaluationDate)
    # TODO: add a check for a remainder
    ProjectionInterval = new_interval(x$EvaluationDate, AsOfDate)
    ProjectionIntervals = suppressMessages(ProjectionInterval / DevelopmentInterval)
    aList = replicate(ProjectionIntervals, x, simplify=FALSE)
    x = do.call("rbind", aList)
    DevIntervals = (1:ProjectionIntervals) * DevelopmentInterval
    x$EvaluationDate = x$EvaluationDate + DevIntervals
    # TODO: zero out stochastic columns
    x
  })
  
  df = do.call("rbind", mojo)
  
  df$DevelopmentLag = as.period(new_interval(df$OriginPeriodStart, df$EvaluationDate + days(1)))
  df$DevInteger = df$DevelopmentLag / DevelopmentInterval
  
  priors = GetPriorNames(objTriangle@StochasticMeasures)
  cumuls = GetCumulativeNames(objTriangle@StochasticMeasures)
  df[, priors] = df[, cumuls]
  
  df
}