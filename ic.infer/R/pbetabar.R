pbetabar <- function (x, df1, df2, wt) 
{
    if (x <= 0) {
        return(0)
    }
    if (df2==0) stop("Error: degrees of freedom (shape2) must not be 0.")
    zed <- df1 == 0
    cdf <- ifelse(any(zed), wt[zed], 0)
    cdf <- cdf + sum(pbeta(x, df1[!zed], df2) * wt[!zed])
    return(cdf)
}