`compcol` <-
function (t, y, params) 
{
    p1 <- y[1]
    p2 <- y[2]
    with(as.list(params), {
        dp1.dt <- c1 * p1 * (1 - p1) - m1 * p1
        dp2.dt <- c2 * p2 * (1 - p1 - p2) - m2 * p2 - c1 * p1 * 
            p2
        return(list(c(dp1.dt, dp2.dt)))
    })
}
