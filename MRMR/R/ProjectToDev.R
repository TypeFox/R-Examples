#' @include TriangleModel.R
#' 

ProjectToDev = function(objTriangleModel, lOriginYears, MaxDev)
{
  objTriangle = objTriangleModel@Triangle
  DevelopmentInterval = objTriangle@DevelopmentInterval
  
  mojo = lapply(lOriginYears, function(x){
    ProjectionIntervals = MaxDev - x$DevInteger
    aList = replicate(ProjectionIntervals, x, simplify=FALSE)
    x = do.call("rbind", aList)
    DevIntervals = (1:ProjectionIntervals) * DevelopmentInterval
    x$EvaluationDate = x$EvaluationDate + DevIntervals
    # TODO: zero out stochastic columns
    x
  })
  
  df = do.call("rbind", mojo)
  
  df$DevelopmentLag = as.period(new_interval(df$OriginPeriodStart, df$EvaluationDate + days(1)))
  df$DevInteger = round(df$DevelopmentLag / DevelopmentInterval, 0)
  
  priors = GetPriorNames(objTriangle@StochasticMeasures)
  cumuls = GetCumulativeNames(objTriangle@StochasticMeasures)
  df[, priors] = df[, cumuls]
  
  df
}