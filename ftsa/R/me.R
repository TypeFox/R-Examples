me = function(forecast, true)
{
   if (length(forecast) != length(true))
       stop("ME: the lengths of input vectors must be the same.")
   err = mean(true - forecast)
   return(round(err, 6))
}
