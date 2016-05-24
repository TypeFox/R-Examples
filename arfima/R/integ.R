# Integrates a series by dint and dseas.

integ <- function(z, zinit = NULL, dint = 0, dseas = 0, period = 0) {
    n <- length(z)
    idcap <- dint + dseas * period
    if (idcap == 0) 
        return(z)
    ncap <- n + idcap
    if (length(zinit) != idcap) 
        stop("zinit not correct length wrt differencing parameters")
    fun <- function(z, ncap, n, dint, dseas, period, zinit, idcap) .Fortran("INTEGD", z = as.double(c(z, 
        rep(0, idcap))), ncap = as.integer(ncap), n = as.integer(n), dint = as.integer(dint), 
        dseas = as.integer(dseas), period = as.integer(period), zinit = as.double(zinit), 
        idcap = as.integer(idcap))
    return(fun(z, ncap, n, dint, dseas, period, zinit, idcap)$z)
} 
