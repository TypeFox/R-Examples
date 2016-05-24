error.den <- function(band, eps, res.data)
{
    res = as.numeric(res.data)
    data.num = length(res)
    std.res = sd(res)
    epsilon = (res - mean(res))/std.res
    eps.std = (eps - mean(res))/std.res
    
    tem  = (eps.std - epsilon)/band
    tem2 = dnorm(tem)/band
    
    hatf = (sum(tem2)/data.num)/std.res
    return(hatf)
}
