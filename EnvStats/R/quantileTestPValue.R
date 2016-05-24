quantileTestPValue <-
function (m, n, r, k, exact.p = TRUE) 
{
    if (!is.vector(m, mode = "numeric") || !is.vector(n, mode = "numeric") || 
        !is.vector(r, mode = "numeric") || !is.vector(k, mode = "numeric")) 
        stop("'m', 'n', 'r',  and 'k' must be numeric vectors")
    if (any(m != trunc(m)) || any(n != trunc(n)) || any(r != 
        trunc(r)) || any(k != trunc(k)) || any(m < 1) || any(n < 
        1) || any(r < 1) || any(k < 1)) 
        stop("All values of 'm', 'n', 'r', and 'k' must be positive integers")
    arg.mat <- cbind.no.warn(m = as.vector(m), n = as.vector(n), 
        r = as.vector(r), k = as.vector(k))
    for (i in c("m", "n", "r", "k")) assign(i, arg.mat[, i])
    N <- m + n
    if (any(r < 2) || any(r > N)) 
        stop("All values of 'r' must be greater than 1 and less than or equal to 'm+n'")
    if (any(k < 0) || any(k > r)) 
        stop(paste("All values of 'k' must be greater than or equal to 0", 
            "and less than or equal to 'r'"))
    if (exact.p) 
        p.val <- 1 - phyper(q = k - 1, m = m, n = n, k = r)
    else p.val <- 1 - pnorm((k - (m * r)/N - 0.5)/sqrt((m * n * 
        r * (N - r))/(N^2 * (N - 1))))
    names(p.val) <- NULL
    p.val
}
