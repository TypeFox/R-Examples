`mdse` <- function(forecast, true)
{
    if (length(forecast) != length(true))
        stop("MSE: the lengths of input vectors must be the same.")
    err = median((true - forecast)^2)
    return(round(err, 4))
}

