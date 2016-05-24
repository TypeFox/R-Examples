`mdape` <- function(forecast, true)
{
     if (length(forecast) != length(true))
         stop("MDAPE: the lengths of input vectors must be the same.")
     err = median(100 * abs((true - forecast) / true))
     return(round(err, 4))
}
