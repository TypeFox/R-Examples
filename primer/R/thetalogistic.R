`thetalogistic` <-
function (times, y, parms) 
{
    n <- y[1]
    with(as.list(parms), {
        dN.dt <- r * n * (1 - (alpha * n)^theta)
        return(list(c(dN.dt)))
    })
}
