qemq <-
function (q, ncp = 0, s = 1) 
{
    if (s <= 0) 
        stop("scaling parameter must be strictly >0.")
    (2 * q - 1)/sqrt(2 * q * (1 - q)) * s + ncp
}
