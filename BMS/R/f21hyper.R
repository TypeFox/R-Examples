f21hyper <-
function (a, b, c, z) 
{
    if ((length(a) != 1) | (length(b) != 1) | (length(c) != 1) | 
        (length(z) != 1)) 
        stop("All function arguments need to be scalars")
    if ((a < 0) | (b < 0) | (c < 0)) 
        stop("Arguments a, b, and c need to be non-negative")
    if ((z > 1) | (z <= (-1))) 
        stop("Argument z needs to be between -1 and 1")
    nmx = max(100, 3 * floor(((a + b) * z - c - 1)/(1 - z)))
    if (nmx > 10000) 
        warning("Power series probably does not converge")
    serie = 0:nmx
    return(1 + sum(cumprod((a + serie)/(c + serie) * (b + serie)/(1 + 
        serie) * z)))
}
