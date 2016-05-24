aovPower <-
function (n.vec, mu.vec = rep(0, length(n.vec)), sigma = 1, alpha = 0.05) 
{
    if (!is.vector(n.vec, mode = "numeric") || !is.vector(mu.vec, 
        mode = "numeric") || !is.vector(sigma, mode = "numeric") || 
        !is.vector(alpha, mode = "numeric")) 
        stop("'n.vec', 'mu.vec', 'sigma', and 'alpha' must be numeric vectors.")
    if (!all(is.finite(n.vec)) || !all(is.finite(mu.vec)) || 
        !all(is.finite(sigma)) || !all(is.finite(alpha))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.vec', 'mu.vec', 'sigma', or 'alpha'"))
    if ((n.grps <- length(n.vec)) < 2 || any(n.vec < 2)) 
        stop(paste("'n.vec' must have at least 2 elements,", 
            "and all values of 'n.vec' must be greater than or equal to 2."))
    if (length(mu.vec) != n.grps) 
        stop("'mu.vec' must be the same length as 'n.vec'")
    if (any(sigma <= 0)) 
        stop("All values of 'sigma' must be positive.")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be between 0 and 1.")
    arg.mat <- cbind.no.warn(sigma = as.vector(sigma), alpha = as.vector(alpha))
    for (i in c("sigma", "alpha")) assign(i, arg.mat[, i])
    N <- length(sigma)
    power.vec <- numeric(N)
    for (i in 1:N) power.vec[i] <- aovPower.scalar(n.vec = n.vec, 
        mu.vec = mu.vec, sigma = sigma[i], alpha = alpha[i])
    names(power.vec) <- NULL
    power.vec
}
