#' PlotModelFactors
#' 
#' @description
#' This function will plot the model factors associated with a triangle model.
#' 
#' @param
#' objTriangleModel A TriangleModel object
#' 
#' @return
#' A ggplot2 plot object
#' 
#' @seealso \code{\link{PlotModelFactors}}
#' 
#' @export PlotModelFactors
#' @include TriangleModel.R
#' 
PlotModelFactors = function(objTriangleModel)
{
  
  fit = objTriangleModel@Fit
  dfCoef = as.data.frame(summary(fit)$coefficients)
  colnames(dfCoef)
  
  x = GetX(objTriangleModel)
  
  dfY = apply(dfCoef, 1, function(y){
    ret = dnorm(x, mean = y[1], sd = y[2])
  })
  
  dfY = as.data.frame(dfY)
  colnames(dfY) = paste0("b", 1:ncol(dfY))
  mdf = suppressMessages(melt(dfY, value.name = "y"))
  mdf$x = rep(x, ncol(dfY))
  
  y = NULL
  variable = NULL
  plt = ggplot(mdf, aes(x, y, col=variable)) + geom_line()
  
  plt
  
}

GetX = function(objTriangleModel)
{
  fit = objTriangleModel@Fit
  
  xlow = min(confint(fit))
  xhigh = max(confint(fit))
  
  x = seq(xlow, xhigh, length.out = 200)
  
  x
}