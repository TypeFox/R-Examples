linearTrendTestN <-
function (slope.over.sigma, alpha = 0.05, power = 0.95, alternative = "two.sided", 
    approx = FALSE, round.up = TRUE, n.max = 5000, tol = 1e-07, 
    maxiter = 1000) 
{
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(slope.over.sigma, mode = "numeric") || !is.vector(alpha, 
        mode = "numeric") || !is.vector(power, mode = "numeric")) 
        stop("'slope.over.sigma', 'alpha', and 'power' must be numeric vectors.")
    if (!all(is.finite(slope.over.sigma)) || !all(is.finite(alpha)) || 
        !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'slope.over.sigma', 'alpha', or 'power'"))
    if (any(abs(slope.over.sigma) < .Machine$double.eps)) 
        stop("All values of 'slope.over.sigma' must be non-zero")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power <= alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than or equal to", 
            "the corresponding elements of 'alpha' and less than 1"))
    if (alternative == "greater" && any(slope.over.sigma <= 0)) 
        stop("When alternative='greater', all values of 'slope.over.sigma' must be positive.")
    if (alternative == "less" && any(slope.over.sigma >= 0)) 
        stop("When alternative='less', all values of 'slope.over.sigma' must be negative.")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        3) 
        stop("'n.max' must be a positive integer greater than 2")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    arg.mat <- cbind.no.warn(power = as.vector(power), slope.over.sigma = as.vector(slope.over.sigma), 
        alpha = as.vector(alpha))
    N <- nrow(arg.mat)
    n.vec <- numeric(N)
    for (i in c("power", "slope.over.sigma", "alpha")) assign(i, 
        arg.mat[, i])
    alt.fac <- ifelse(alternative == "two.sided", 2, 1)
    fcn.for.root <- function(n, power, slope.over.sigma, alpha, 
        alternative, approx) {
        power - linearTrendTestPower(n = n, slope.over.sigma = slope.over.sigma, 
            alpha = alpha, alternative = alternative, approx = approx)
    }
    power.3 <- linearTrendTestPower(n = 3, slope.over.sigma = slope.over.sigma, 
        alpha = alpha, alternative = alternative, approx = approx)
    power.n.max <- linearTrendTestPower(n = n.max, slope.over.sigma = slope.over.sigma, 
        alpha = alpha, alternative = alternative, approx = approx)
    for (i in 1:N) {
        power.i <- power[i]
        power.3.i <- power.3[i]
        if (power.3.i >= power.i) 
            n.vec[i] <- 3
        else {
            power.n.max.i <- power.n.max[i]
            if (power.n.max.i < power.i) {
                n.vec[i] <- NA
                warning("Error in search algorithm.  Try increasing the value of the argument 'n.max'")
            }
            else {
                slope.over.sigma.i <- slope.over.sigma[i]
                alpha.i <- alpha[i]
                n.vec[i] <- uniroot(fcn.for.root, lower = 3, 
                  upper = n.max, f.lower = power.i - power.3.i, 
                  f.upper = power.i - power.n.max.i, power = power.i, 
                  slope.over.sigma = slope.over.sigma.i, alpha = alpha.i, 
                  alternative = alternative, approx = approx, 
                  tol = tol, maxiter = maxiter)$root
            }
        }
    }
    if (round.up) 
        n.vec <- ceiling(n.vec)
    n.vec
}
