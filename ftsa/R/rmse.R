`rmse` <- function(forecast, true)
{
    if (length(forecast) != length(true))
        stop("RMSE: the lengths of input vectors must be the same.")
    err = sqrt(mean((true - forecast)^2))
    return(round(err, 6))
}

