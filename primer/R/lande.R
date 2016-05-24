`lande` <-
function (t, y, parms) 
{
    p <- y[1]
    with(as.list(parms), {
        dp <- ci * p * (1 - D - p) - e * p
        return(list(dp))
    })
}
