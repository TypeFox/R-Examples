mrae = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("MRAE: the lengths of input vectors must be the same.")
  ferror = (true - forecast)
  ferrorbench = (true - forecastbench)
  relativerror = mean(abs(ferror / ferrorbench))
  return(round(relativerror, 6))
}
