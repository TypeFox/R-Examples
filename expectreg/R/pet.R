pet <-
function (e, df) 
{
    u = (df + e^2) * dt(e, df)/(1 - df) - e * pt(e, df)
    asy = u/(2 * u + e)
    return(asy)
}
