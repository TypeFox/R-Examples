`smape` <-
function(forecast, true)
{
     if (length(forecast) != length(true))
         stop("SMAPE: the lengths of input vectors must be the same.")
     err = mean(200 * (true - forecast) / (true + forecast))
     return(round(err, 6))
}
