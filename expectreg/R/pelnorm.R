pelnorm <-
function (e, meanlog = 0, sdlog = 1) 
{
    u = exp(meanlog + 0.5 * sdlog^2) * (1 - pnorm((meanlog + 
        sdlog^2 - log(e))/sdlog)) - e * plnorm(e, meanlog, sdlog)
    asy = u/(2 * u + e - exp(meanlog + 0.5 * sdlog^2))
    return(asy)
}
