`SIRbd` <-
function (t, y, p) 
{
    S <- y[1]
    I <- y[2]
    R <- y[3]
    with(as.list(p), {
        dS.dt <- b * (S + I + R) - B * I * S - m * S
        dI.dt <- B * I * S - g * I - m * I
        dR.dt <- g * I - m * R
        return(list(c(dS.dt, dI.dt, dR.dt)))
    })
}
