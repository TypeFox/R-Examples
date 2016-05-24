`mape` <- function(forecast, true)
{
     if (length(forecast) != length(true))
         stop("MAPE: the lengths of input vectors must be the same.")
     err = mean(100 * abs((true - forecast) / true))
     return(round(err, 6))
}

