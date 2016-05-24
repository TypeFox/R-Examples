ponderer.classe <-
function(data)
{
    n <- nrow(data)
    if (n == 1)
    {
        return(exp(1)/(1+exp(1)))
    }

    vx <- var(data$x)*(n-1)/n
    vy <- var(data$y)*(n-1)/n

    cxy <- cov(data$x, data$y)*(n-1)/n

    delta <- vx*vy -cxy**2
    pond <- exp(n/(n+delta))/(1+exp(n/(n+delta)))

    return(pond)
}
