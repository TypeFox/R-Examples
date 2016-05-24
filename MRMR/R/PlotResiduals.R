#' @title PlotResiduals
#' 
#' @export
#' 
#' @include TriangleModel.R
#' 
#' @description
#' This will produce a 2x2 set of residual graphs.
#' 
#' @param objTriangleModel A TriangleModel object
#' 
#' @details
#' This function will produce four charts. 
#' 
#' @return 
#' This function does not return a value.
#' 
PlotResiduals = function(objTriangleModel){
  fit = objTriangleModel@Fit

  residuals = rstandard(fit)
  fittedValues = fitted(fit)
  
  predictor = objTriangleModel@Predictor
  whichRows = !is.na(objTriangleModel@ModelData[predictor])
  dfModelData = objTriangleModel@ModelData[whichRows,]
  
  devLag = dfModelData$DevInteger
  originPeriod = dfModelData$OriginPeriodStart
  calendarPeriod = dfModelData$CalendarPeriodStart
  dfModelData$residuals = residuals

#  op = par(mfrow=c(2,2), oma=c(0,0,2,0))
  op = par(mfrow=c(2,2))
  
  ylow = min(residuals, -3)
  ymax = max(residuals, 3)
#   ylow = -3
#   ymax = 3
  
  plot(x = fittedValues, y = residuals, col="blue", pch=16, xlab="Fitted values", ylab="Std. Resid.", xaxt="n", ylim=c(ylow, ymax))
  xTicks = pretty(fittedValues, 8)
  axis(side=1, at=xTicks, labels=format(xTicks, scientific=FALSE, big.mark=","))
  abline(0, 0)
  abline(3, 0, lty = "dashed")
  abline(-3, 0, lty = "dashed")

  plot(devLag, residuals, col="blue", pch=16, xlab="Development Interval", ylab="Std. Resid.", ylim=c(ylow, ymax))
  abline(0,0)
  abline(3, 0, lty = "dashed")
  abline(-3, 0, lty = "dashed")
  factors = factor(devLag)
  means = tapply(residuals, factors, mean)
  lines(levels(factors), means, col="red", type="l")

  plot(originPeriod, residuals, col="blue", pch=16, xlab="Origin Period", ylab="Std. Resid.", ylim=c(ylow, ymax))
  abline(0,0)
  abline(3, 0, lty = "dashed")
  abline(-3, 0, lty = "dashed")
  OriginPeriodStart = NULL
  means = daply(dfModelData, .(OriginPeriodStart), function(df) mean(df$residuals))
  originPeriods = unique(originPeriod)
  lines(originPeriods, means, col="red", type="l")

  plot(calendarPeriod, residuals, col="blue", pch=16, xlab="Calendar Period", ylab="Std. Resid.", ylim=c(ylow, ymax))
  abline(0,0)
  abline(3, 0, lty = "dashed")
  abline(-3, 0, lty = "dashed")
  CalendarPeriodStart = NULL
  means = daply(dfModelData, .(CalendarPeriodStart), function(df) mean(df$residuals))
  calPeriods = unique(calendarPeriod)
  lines(calPeriods, means, col="red", type="l")
   
#  mtext(Title, side=3, outer=TRUE, line=0, adj=0)

  par(op)
  
}