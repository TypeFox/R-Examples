peunif <-
function (e, min = 0, max = 1) 
{
    u = (e^2 - min^2)/2/(max - min) - e * punif(e, min = min, 
        max = max)
    asy = u/(2 * u + e - (min + max)/2)
    return(asy)
}
