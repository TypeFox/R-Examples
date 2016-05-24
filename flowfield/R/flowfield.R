flowfield <- function (t, y, steps,plot) {
  #*****************************************************************************
  #  Flow field forecasting draws information from an interpolated flow field of 
  #  the observed time series to incrementally build a forecast. The time series 
  #  need not have uniformly spaced observations. Flow field forecasting works 
  #  best on relatively long time series (i.e. > 1000 observations) where 
  #  forecasts must be made autonomously. 
  #
  #   Input: t - time series observation times 
  #          y - time series response values
  #          steps - Number of steps (1-10) to forecast. Forecasts occur in knot 
  #                  intervals of the penalized spline regression. Knots are 
  #                  evenly spaced within the range of data appoximately one knot 
  #                  for every 10 data points
  #          plot -  If a plot is required, set plot = TRUE
  #
  #   Output:  Forecast values, prediction errors at each step, and plot if 
  #            selected
  #******************************************************************************
  
  skeleton <- psr(t,y)
  if ((steps > 10) || (steps < 0)) {
    stop("Warning:  Steps must be between 0 and 10", file="")
  }
  else {fcast <- forecast(skeleton,steps)
        if (plot==TRUE) {ffplot(t,y,skeleton,fcast$forecast,fcast$error)}
        output <- write.table(fcast,file="",row.names=FALSE,quote=FALSE,sep = '\t\t')
  }
  return(output)
}