remq <-
function (n, ncp = 0, s = 1) 
{
    if (s <= 0) 
        stop("scaling parameter must be strictly >0.")
    x = runif(n)
    (2 * x - 1)/sqrt(2 * x * (1 - x)) * s + ncp
}
