peemq <-
function (e, ncp = 0, s = 1) 
{
    z = (e - ncp)/s
    u = 0 * e
    for (k in 1:length(e)) u[k] = integrate(function(x) x * demq(x), 
        -Inf, z[k])$value - z[k] * pemq(z[k])
    asy = u/(2 * u + z)
    return(asy)
}
