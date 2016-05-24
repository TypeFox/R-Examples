pebeta <-
function (e, a = 1, b = 1) 
{
    u = gamma(1 + a) * gamma(a + b)/gamma(a)/gamma(1 + a + b) * 
        pbeta(e, 1 + a, b) - e * pbeta(e, a, b)
    asy = u/(2 * u + e - a/(a + b))
    return(asy)
}
