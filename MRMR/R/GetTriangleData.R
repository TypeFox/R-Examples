#' @title GetTriangleData
#' @export GetTriangleData
#' @include Triangle.R
#' 
#' @description
#' This function will return data values from a triangle.
#' 
#' @param Triangle A Triangle object
#' @param OriginPeriodStart A vector of origin years. This parameter may be null.
#' @param DevInteger A vector of development integers. This parameter may be null.
#' @param EvaluationDate A vector of evaluation dates. This parameter may be null.
#' @param Measure A character vector with the names of measures to return.
#' 
#' @return A data frame
#' 
GetTriangleData = function(Triangle, OriginPeriodStart = NULL, DevInteger = NULL, EvaluationDate = NULL, Measure){
  df = Triangle@TriangleData
  
  if (!is.null(OriginPeriodStart)){
    tz(OriginPeriodStart) = tz(df$OriginPeriodStart)
    df = df[df$OriginPeriodStart %in% OriginPeriodStart, ]
  } 
  
  if (!is.null(DevInteger)){
    df = df[df$DevInteger %in% DevInteger, ]
  } 
  
  if (!is.null(EvaluationDate)){
    tz(EvaluationDate) = tz(df$EvaluationDate)
    df = df[df$EvaluationDate %in% EvaluationDate, ]
  } 
  
  df = df[, colnames(df) %in% Measure]
  
  row.names(df) = NULL
  
  df
}

mojo = TestDataFrame()
myTriangle = newTriangle(TriangleData = mojo
                         , OriginPeriods = AccidentYear
                         , DevelopmentLags = Month
                         , Cumulative = TRUE
                         , StochasticMeasures = c("Paid")
                         , StaticMeasures = c("EP")
                         , Verbose = FALSE)
cazart = GetTriangleData(myTriangle, OriginPeriodStart = mdy("1-1-2002"), Measure = c("IncrementalPaid", "PriorPaid"))
cazart

cazart = GetTriangleData(myTriangle, OriginPeriodStart = mdy("1-1-2002"), DevInteger = 2, Measure = c("IncrementalPaid", "PriorPaid"))
cazart

cazart = GetTriangleData(myTriangle, DevInteger = 2, Measure = c("IncrementalPaid", "PriorPaid"))
cazart

cazart = GetTriangleData(myTriangle, EvaluationDate = ymd("2004-12-31"), Measure = c("IncrementalPaid", "PriorPaid"))
cazart
