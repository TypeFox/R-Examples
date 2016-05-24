pchibar <- function (x, df, wt) 
{
    ## acknowlegment: taken from package ibdreg
    if (x <= 0) {
        return(0)
    }
    zed <- df == 0
    cdf <- ifelse(any(zed), wt[zed], 0)
    cdf <- cdf + sum(pchisq(x, df[!zed]) * wt[!zed])
    return(cdf)
}
