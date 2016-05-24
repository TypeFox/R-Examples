`predpreyRM` <-
function (t, y, p) 
{
    H <- y[1]
    P <- y[2]
    with(as.list(p), {
        dH.dt <- b * H * (1 - alpha * H) - w * P * H/(D + H)
        dP.dt <- e * w * P * H/(D + H) - s * P
        return(list(c(dH.dt, dP.dt)))
    })
}
