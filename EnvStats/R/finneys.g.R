finneys.g <-
function (m, z, n.terms.inc = 10, iter = 100, tol = .Machine$double.eps) 
{
    if (any(!is.finite(m)) || any(!is.finite(z))) 
        stop(paste("Missing (NA), infinite (-Inf, Inf), and undefined (NaN) values", 
            "not allowed in 'm' or 'z'"))
    if (any(m < 1) || any(m != trunc(m))) 
        stop("All values of 'm' must be integers greater than or equal to 1")
    if (any(n.terms.inc < 3) || any(n.terms.inc != trunc(n.terms.inc))) 
        stop("All values of 'n.terms.inc' must be integers greater than or equal to 3")
    if (any(iter < 1) || any(iter != trunc(iter))) 
        stop("All values of 'iter' must be integers greater than or equal to 1")
    if (any(tol < .Machine$double.eps)) 
        stop("All values of 'tol' must be at least as large as .Machine$double.eps")
    arg.mat <- cbind.no.warn(m = as.vector(m), z = as.vector(z), 
        n.terms.inc = as.vector(n.terms.inc), iter = as.vector(iter), 
        tol = as.vector(tol))
    for (i in c("m", "z", "n.terms.inc", "iter", "tol")) assign(i, 
        arg.mat[, i])
    n <- length(m)
    vec <- numeric(n)
    for (i in 1:n) vec[i] <- finneys.g.scalar(m[i], z[i], n.terms.inc[i], 
        iter[i], tol[i])
    vec
}
