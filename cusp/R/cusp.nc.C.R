`cusp.nc.C` <-
function (alpha, beta, subdivisions = 100, rel.tol = .Machine$double.eps^0.25, 
    abs.tol = rel.tol, stop.on.error = TRUE, aux = NULL, keep.order = TRUE) 
{
    limit <- as.integer(subdivisions)
    if (limit < 1 || (abs.tol <= 0 && rel.tol < max(50 * .Machine$double.eps, 
        5e-29))) 
        stop("invalid parameter values")
    lower <- -Inf
    upper <- Inf
    inf <- 2
    bound <- 0
    wk <- .External("cuspnc", as.double(alpha), as.double(beta), 
        as.double(bound), as.integer(inf), as.double(abs.tol), 
        as.double(rel.tol), limit = limit)
    wk
}

