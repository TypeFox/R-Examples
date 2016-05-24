`mpe` <- function(forecast, true)
{
     if (length(forecast) != length(true))
         stop("MPE: the lengths of input vectors must be the same.")
     err = mean(100 * ((true - forecast) / true))
     return(round(err, 6))
}

