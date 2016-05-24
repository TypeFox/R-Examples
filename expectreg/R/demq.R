demq <-
function (z, ncp = 0, s = 1) 
{
    if (s <= 0) 
        stop("scaling parameter must be strictly > 0.")
    x = (z - ncp)/s
    (1/s)/(2 + x^2)^1.5
}
