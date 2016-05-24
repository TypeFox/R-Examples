pemq <-
function (z, ncp = 0, s = 1) 
{
    if (s <= 0) 
        stop("scaling parameter must be strictly >0.")
    x = (z - ncp)/s
    0.5 + 0.5 * sign(x) * sqrt(1 - 2/(2 + x^2))
}
