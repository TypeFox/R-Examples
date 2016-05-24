`rmdspe` <- function(forecast, true)
{
     if (length(forecast) != length(true))
         stop("RMDSPE: the lengths of input vectors must be the same.")
     err = sqrt(median((100 * (true - forecast) / true)^2))
     return(round(err, 6))
}
