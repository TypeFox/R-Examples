# Summary Function to calclulate mean and quantiles
# for genrated forecasts

"summary.forecast" <- function(object, probs = c(0.16, 0.84), ...)
#identify the type of forecast
{
  if(inherits(object, "forecast.VAR")){
    return(summary.forecast.VAR(object, probs = probs))
    }

  if(inherits(object, "forecast.BVAR")){
    return(summary.forecast.BVAR(object, probs = probs))
    }

  if(inherits(object, "forecast.BSVAR")){
    return(summary.forecast.BSVAR(object, probs = probs))
    }

}
"summary.forecast.VAR" <- function(object, probs, ...)
{
  #Mean
  mean.forecast.Var <- apply(object, 2, mean)
  #Quantile
  if (missing(probs)){
  quantile.forecast.Var <- apply(object, 2, quantile, probs = c(0.16, 0.84))
  }
  else{
  quantile.forecast.Var <- apply(object, 2, quantile, probs)
  }
  cat("Summary\nMean\n")
  print(mean.forecast.Var)
  cat("\nQuantiles\n")
  print(quantile.forecast.Var)
}

"summary.forecast.BVAR" <- function(object, probs, ...)
{
  output <- summary.forecast.VAR(object, probs)
}

"summary.forecast.BSVAR" <- function(object, probs, ...)
{
  output <- summary.forecast.VAR(object, probs)
}
