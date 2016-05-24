aovN <-
function (mu.vec, sigma = 1, alpha = 0.05, power = 0.95, round.up = TRUE, 
    n.max = 5000, tol = 1e-07, maxiter = 1000) 
{
    if (!is.vector(mu.vec, mode = "numeric") || !is.vector(sigma, 
        mode = "numeric") || !is.vector(alpha, mode = "numeric") || 
        !is.vector(power, mode = "numeric")) 
        stop(paste("'mu.vec', 'sigma, 'alpha', and 'power'", 
            "must be numeric vectors."))
    if (!all(is.finite(mu.vec)) || !all(is.finite(sigma)) || 
        !all(is.finite(alpha)) || !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'mu.vec', 'sigma', 'alpha', or 'power'"))
    if ((n.grps <- length(mu.vec)) < 2 || length(unique(mu.vec)) == 
        1) 
        stop("'mu.vec' must have at least 2 distinct values")
    if (any(sigma <= 0)) 
        stop("All values of 'sigma' must be positive.")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power <= alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than the", 
            "corresponding elements of 'alpha' and less than 1"))
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    arg.mat <- cbind.no.warn(sigma = as.vector(sigma), alpha = as.vector(alpha), 
        power = as.vector(power))
    for (i in c("sigma", "alpha", "power")) assign(i, arg.mat[, 
        i])
    N <- length(sigma)
    n.vec <- numeric(N)
    for (i in 1:N) n.vec[i] <- aovN.scalar(mu.vec = mu.vec, sigma = sigma[i], 
        alpha = alpha[i], power = power[i], round.up = round.up, 
        n.max = n.max, tol = tol, maxiter = maxiter)
    names(n.vec) <- NULL
    n.vec
}
