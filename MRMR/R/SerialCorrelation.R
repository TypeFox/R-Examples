GetPriorResiduals = function(df, Axis1, Axis2){
  rowOrder = order(df[Axis1], df[Axis2])
  df = df[rowOrder, ]
  df$PriorResidual = c(NA, df$Residual[1:(nrow(df)-1)])
  
  # This is a hack to get around the fact that there is no min function in lubridate
  # We order the data frame by Axis2 and NA everything equal to the first row
  checkVal = df[order(df[Axis2]), Axis2]
  naRow = (df[,Axis2] == checkVal[1])
  df$PriorResidual[naRow] = NA
  
  df = df[!naRow, ]
  row.names(df) = NULL
  
  df
}

#' Fit the serial correlation in a triangle
#' @export 
#' 
#' @include TriangleModel.R
#' 
#' @param objTriangleModel A Triangle model
FitSerialCorrelation = function(objTriangleModel){
  df = GetPriorResiduals(objTriangleModel@ModelData, "DevInteger", "OriginPeriodStart")
  df$EvalGroup = as.factor(df$EvaluationDate)
  
  fit = lm(Residual ~ 0 + PriorResidual:EvalGroup, data = df)
  
  lstReturn = list(df = df, fit = fit)
}

#' Plot the serial correlation in a triangle
#' @export 
#' 
#' @include TriangleModel.R
#' 
#' @param objTriangleModel A Triangle model
plotSerialCorrelation = function(objTriangleModel){
  lstReturn = FitSerialCorrelation(objTriangleModel)
  df = lstReturn$df
  fit = lstReturn$fit
  PriorResidual = NULL
  Residual = NULL
  EvalGroup = NULL
  plt = ggplot(df, aes(x = PriorResidual, y = Residual, group = EvalGroup, colour = EvalGroup)) + geom_point()
  plt = plt + stat_smooth(method = lm, se = FALSE)
  print(plt)
  
  lstReturn
}

