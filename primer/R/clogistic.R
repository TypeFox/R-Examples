`clogistic` <-
function (times, y, parms) 
{
    n <- y[1]
    r <- parms[1]
    alpha <- parms[2]
    dN.dt <- r * n * (1 - alpha * n)
    return(list(c(dN.dt)))
}
