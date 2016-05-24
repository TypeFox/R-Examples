sse = function (forecast, true) 
{
    if (length(forecast) != length(true)) 
        stop("SSE: the lengths of input vectors must be the same.")
    err = sum((true - forecast)^2)
    return(round(err, 6))
}