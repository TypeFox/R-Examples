`dcusp` <-
function (y, alpha, beta) 
{
    f <- function(x, a, b) {
        x2 = x^2
        exp(a * x + b * x2/2 - x2^2/4)
    }
    C <- if (is.loaded("cuspnc")) 
        cusp.nc.c(alpha, beta)
    else cusp.nc(alpha, beta)
    f(y, alpha, beta)/C
}

