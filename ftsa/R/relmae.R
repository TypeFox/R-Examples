relmae = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("RelMAE: the lengths of input vectors must be the same.")
  ferror = mean(abs(true - forecast))
  ferrorbench = mean(abs(true - forecastbench))  
  relativerror = ferror / ferrorbench
  return(round(relativerror, 6))
}
