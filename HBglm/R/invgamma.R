dinvgamma <- function(y, shape, scale, log=FALSE)
{
    x     <- 1 / y
    scale <- 1 / scale
    s     <- dgamma(x, shape, scale=scale, log=log)
    if (log) return(s + 2 * log(x))
    return(s * (x^2))
}

rinvgamma <- function(n, shape, scale)
{
    return(1/rgamma(n, shape, rate=scale))
}
