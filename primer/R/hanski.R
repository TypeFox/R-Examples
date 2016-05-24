`hanski` <-
function (t, y, parms) 
{
    p <- y[1]
    with(as.list(parms), {
        dp <- ci * p * (1 - p) - e * p * (1 - p)
        return(list(dp))
    })
}
