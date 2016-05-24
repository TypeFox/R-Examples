`SIR` <-
function (t, y, p) 
{
    {
        S <- y[1]
        I <- y[2]
        R <- y[3]
    }
    with(as.list(p), {
        dS.dt <- -B * I * S
        dI.dt <- B * I * S - g * I
        dR.dt <- g * I
        return(list(c(dS.dt, dI.dt, dR.dt)))
    })
}
