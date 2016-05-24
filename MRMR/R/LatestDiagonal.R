#' @title LatestDiagonal
#' @export LatestDiagonal
#' @include Triangle.R
#' 
#' @description
#' This function will return all of the values for the most recent evaluation date. Note 
#' that this applies for each origin period individually. For example, if some origin periods have an 
#' evaluation at December 31, 2010, but others only have evaluations at December 31, 2009, the data frame 
#' which is returned will have two different evaluation dates present.
#' 
#' @param
#' x a data frame or a triangle
#' 
#' @return
#' A data frame
#' 
LatestDiagonal = function(x)
{
  if (class(x) == "Triangle") x = x@TriangleData
  
  lOriginYear = dlply(x, "OriginPeriodStart")
  lOriginYear = lapply(lOriginYear, function(y) {
    whichRows = y$EvaluationDate == max(y$EvaluationDate)
    latest = y[whichRows,]
  })
  
  dfLatest = do.call("rbind", lOriginYear)
  
  dfLatest
  
}