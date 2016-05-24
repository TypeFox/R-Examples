relmse = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("RelMAE: the lengths of input vectors must be the same.")
  ferror = mean((true - forecast)^2)
  ferrorbench = mean((true - forecastbench)^2)
  relativerror = ferror / ferrorbench
  return(round(relativerror, 6))
}
