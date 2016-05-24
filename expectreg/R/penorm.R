penorm <-
function (e, m = 0, sd = 1) 
{
    z = (e - m)/sd
    p = pnorm(z)
    d = dnorm(z)
    u = -d - z * p
    asy = u/(2 * u + z)
    return(asy)
}
