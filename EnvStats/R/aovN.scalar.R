aovN.scalar <-
function (mu.vec, sigma = 1, alpha = 0.05, power = 0.95, round.up = TRUE, 
    n.max = 5000, tol = 1e-07, maxiter = 1000) 
{
    if (!is.vector(mu.vec, mode = "numeric")) 
        stop("'mu.vec' must be a numeric vector.")
    if (!is.vector(sigma, mode = "numeric") || length(sigma) != 
        1 || !is.vector(alpha, mode = "numeric") || length(alpha) != 
        1 || !is.vector(power, mode = "numeric") || length(power) != 
        1) 
        stop("'sigma, 'alpha', and 'power' must be numeric scalars.")
    if (!all(is.finite(mu.vec)) || !is.finite(sigma) || !is.finite(alpha) || 
        !is.finite(power)) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'mu.vec', 'sigma', 'alpha', or 'power'"))
    if ((n.grps <- length(mu.vec)) < 2 || length(unique(mu.vec)) == 
        1) 
        stop("'mu.vec' must have at least 2 distinct values")
    if (sigma <= 0) 
        stop("'sigma' must be positive.")
    if (alpha <= 0 || alpha >= 1) 
        stop("'alpha' must be greater than 0 and less than 1")
    if (power <= alpha || power >= 1) 
        stop(paste("All values of 'power' must be greater than the", 
            "corresponding elements of 'alpha' and less than 1"))
    fcn.for.root <- function(n, power, mu.vec, sigma, alpha, 
        n.grps) {
        power - aovPower(n.vec = rep(n, n.grps), mu.vec = mu.vec, 
            sigma = sigma, alpha = alpha)
    }
    power.2 <- aovPower(n.vec = rep(2, n.grps), mu.vec = mu.vec, 
        sigma = sigma, alpha = alpha)
    if (power.2 >= power) 
        ret.val <- 2
    else {
        delta <- max(diff(sort(mu.vec)))
        power.n.max <- aovPower(n.vec = rep(n.max, n.grps), mu.vec = mu.vec, 
            sigma = sigma, alpha = alpha)
        if (power.n.max < power) {
            n <- NA
            warning(paste("Error in search algorithm.", "Try increasing the value of the argument 'n.max'."))
        }
        else {
            n <- uniroot(fcn.for.root, lower = 2, upper = n.max, 
                f.lower = power - power.2, f.upper = power - 
                  power.n.max, power = power, mu.vec = mu.vec, 
                sigma = sigma, alpha = alpha, n.grps = n.grps, 
                tol = tol, maxiter = maxiter)$root
        }
        ret.val <- ifelse(round.up, ceiling(n), n)
    }
    ret.val
}
