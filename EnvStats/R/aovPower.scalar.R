aovPower.scalar <-
function (n.vec, mu.vec = rep(0, length(n.vec)), sigma = 1, alpha = 0.05) 
{
    if (!is.vector(n.vec, mode = "numeric") || !is.vector(mu.vec, 
        mode = "numeric")) 
        stop("'n.vec' and 'mu.vec' must be numeric vectors.")
    if (!is.vector(sigma, mode = "numeric") || length(sigma) != 
        1 || !is.vector(alpha, mode = "numeric") || length(alpha) != 
        1) 
        stop("'sigma' and 'alpha' must be numeric scalars")
    if (!all(is.finite(n.vec)) || !all(is.finite(mu.vec)) || 
        !is.finite(sigma) || !is.finite(alpha)) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.vec', 'mu.vec', 'sigma', or 'alpha'"))
    if ((n.grps <- length(n.vec)) < 2 || any(n.vec < 2)) 
        stop(paste("'n.vec' must have at least 2 elements,", 
            "and all values of 'n.vec' must be greater than or equal to 2."))
    if (length(mu.vec) != n.grps) 
        stop("'mu.vec' must be the same length as 'n.vec'")
    if (sigma <= 0) 
        stop("'sigma' must be positive.")
    if (alpha <= 0 || alpha >= 1) 
        stop("'alpha' must be between 0 and 1.")
    df1 <- n.grps - 1
    df2 <- sum(n.vec) - n.grps
    x <- rep(mu.vec, n.vec)
    A <- factor(rep(1:n.grps, n.vec))
    SS.res <- sum(sapply(split(x, A), function(z) (length(z) - 
        1) * var(z)))
    SS.tot <- (length(x) - 1) * var(x)
    ncp <- (SS.tot - SS.res)/sigma^2
    1 - pf(qf(1 - alpha, df1, df2), df1, df2, ncp = ncp)
}
