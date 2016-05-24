`gotelli` <-
function (t, y, parms) 
{
    p <- y[1]
    with(as.list(parms), {
        dp <- ce * (1 - p) - e * p
        return(list(dp))
    })
}
