`predpreyLV` <-
function (t, y, params) 
{
    H <- y[1]
    P <- y[2]
    with(as.list(params), {
        dH.dt <- b * H - a * P * H
        dP.dt <- e * a * P * H - s * P
        return(list(c(dH.dt, dP.dt)))
    })
}
