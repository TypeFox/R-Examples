`cusp.logLike` <-
function (p, x, verbose = FALSE) 
{
    if (verbose) {
        print(p)
        flush.console()
    }
    z = (x - p[3])/p[4]
    -2 * sum(p[1] * z + p[2] * z^2/2 - z^4/4) + 2 * length(x) * 
        log(p[4] * cusp.nc(p[1], p[2]))
}

