`SIRf` <-
function (t, y, p) 
{
    {
        S <- y[1]
        I <- y[2]
        R <- y[3]
        N <- S + I + R
    }
    with(as.list(p), {
        dS.dt <- -B * I * S/N
        dI.dt <- B * I * S/N - g * I
        dR.dt <- g * I
        return(list(c(dS.dt, dI.dt, dR.dt)))
    })
}
