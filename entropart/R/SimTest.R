as.SimTest <- 
function (RealValue, SimulatedValues) {
  st <- list(RealValue=RealValue, SimulatedValues=SimulatedValues)
  class(st) <- "SimTest"
  return(st)
}


is.SimTest <-
function (x) 
{
  inherits(x, "SimTest")
}


plot.SimTest <- 
function (x, Quantiles = c(0.025, 0.975), ..., colValue = "red", lwdValue = 2, ltyValue = 2, colQuantiles = "black", lwdQuantiles = 1, ltyQuantiles = 2)
{
  plot(stats::density(x$SimulatedValues), ...)
  graphics::abline(v=x$RealValue, col=colValue, lwd=lwdValue, lty=ltyValue)
  for (qt in Quantiles) {
    graphics::abline(v=stats::quantile(x$SimulatedValues, probs = qt), col=colQuantiles, lwd=lwdQuantiles, lty=ltyQuantiles)
  }
}


summary.SimTest <-
function(object, Quantiles = c(0.025, 0.975), ...) 
{ 
  cat("Real value: ", object$RealValue, "\n")
  cat("Quantile in the simulated distribution: ", stats::ecdf(object$SimulatedValues)(object$RealValue), "\n")
  
  cat("Quantiles of simulations:\n")
  for (qt in Quantiles) {
    cat(sprintf("%1.2f%%", 100*qt), ": ", stats::quantile(object$SimulatedValues, probs = qt), "\n")
  }
  cat("Mean simulated value: ", mean(object$SimulatedValues))
  return(invisible(NULL))
}
