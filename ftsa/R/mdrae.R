mdrae = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("MdRAE: the lengths of input vectors must be the same.")
  ferror = (true - forecast)
  ferrorbench = (true - forecastbench)
  relativerror = median(abs(ferror / ferrorbench))
  return(round(relativerror, 6))
}  
