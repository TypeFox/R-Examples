"kinhat" <-
function (pts, lambda, poly, s)
{
    ptsx <- pts[, 1]
    ptsy <- pts[, 2]
    npt <- length(ptsx)
    ns <- length(s)
    s <- sort(s)
    np <- length(poly[, 1])
    polyx <- c(poly[, 1], poly[1, 1])
    polyy <- c(poly[, 2], poly[1, 2])
    hkhat <- rep(0, times = ns)
    klist <- .Fortran("dokinhat", as.double(ptsx), as.double(ptsy),
        as.integer(npt), as.double(lambda), as.double(polyx), as.double(polyy),
        as.integer(np), as.double(s), as.integer(ns), as.double(hkhat), 
        PACKAGE="spatialkernel")
    res <- list(k = as.numeric(klist[[10]]), s = s)
}
