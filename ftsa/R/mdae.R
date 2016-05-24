mdae = function(forecast, true)
{
    if (length(forecast) != length(true))
        stop("MDAE: the lengths of input vectors must be the same.")
    err = median(abs(true - forecast))
    return(round(err, 4))
}

