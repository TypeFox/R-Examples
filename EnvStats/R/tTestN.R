tTestN <-
function (delta.over.sigma, alpha = 0.05, power = 0.95, sample.type = ifelse(!is.null(n2), 
    "two.sample", "one.sample"), alternative = "two.sided", approx = FALSE, 
    n2 = NULL, round.up = TRUE, n.max = 5000, tol = 1e-07, maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(delta.over.sigma, mode = "numeric") || !is.vector(alpha, 
        mode = "numeric") || !is.vector(power, mode = "numeric")) 
        stop("'delta.over.sigma', 'alpha', and 'power' must be numeric vectors.")
    if (!all(is.finite(delta.over.sigma)) || !all(is.finite(alpha)) || 
        !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'delta.over.sigma', 'alpha', or 'power'"))
    if (any(abs(delta.over.sigma) <= 0)) 
        stop("All values of 'delta.over.sigma' must be non-zero")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power <= alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than or equal to", 
            "the corresponding elements of 'alpha' and less than 1"))
    if (alternative == "greater" && any(delta.over.sigma <= 0)) 
        stop("When alternative='greater', all values of 'delta.over.sigma' must be positive.")
    if (alternative == "less" && any(delta.over.sigma >= 0)) 
        stop("When alternative='less', all values of 'delta.over.sigma' must be negative.")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    if (n2.constrained <- sample.type == "two.sample" && !is.null(n2)) {
        if (!is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (!all(is.finite(n2))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in 'n2'"))
        if (any(n2 != trunc(n2)) || any(n2 < 2)) 
            stop("All values of 'n2' must be positive integers larger than 1")
        arg.mat <- cbind.no.warn(power = as.vector(power), delta.over.sigma = as.vector(delta.over.sigma), 
            alpha = as.vector(alpha), n2 = as.vector(n2))
        N <- nrow(arg.mat)
        n.vec <- numeric(N)
        for (i in c("power", "delta.over.sigma", "alpha", "n2")) assign(i, 
            arg.mat[, i])
    }
    else {
        arg.mat <- cbind.no.warn(power = as.vector(power), delta.over.sigma = as.vector(delta.over.sigma), 
            alpha = as.vector(alpha))
        N <- nrow(arg.mat)
        n.vec <- numeric(N)
        for (i in c("power", "delta.over.sigma", "alpha")) assign(i, 
            arg.mat[, i])
    }
    type.fac <- ifelse(sample.type == "two.sample", 2, 1)
    alt.fac <- ifelse(alternative == "two.sided", 2, 1)
    sod2 <- 1/delta.over.sigma^2
    fcn.for.root <- function(n, power, delta.over.sigma, alpha, 
        sample.type, alternative, approx) {
        power - tTestPower(n.or.n1 = n, delta.over.sigma = delta.over.sigma, 
            alpha = alpha, sample.type = sample.type, alternative = alternative, 
            approx = approx)
    }
    power.2 <- tTestPower(n.or.n1 = 2, delta.over.sigma = delta.over.sigma, 
        alpha = alpha, sample.type = sample.type, alternative = alternative, 
        approx = approx)
    power.n.max <- tTestPower(n.or.n1 = n.max, delta.over.sigma = delta.over.sigma, 
        alpha = alpha, sample.type = sample.type, alternative = alternative, 
        approx = approx)
    for (i in 1:N) {
        power.i <- power[i]
        power.2.i <- power.2[i]
        if (power.2.i >= power.i) 
            n.vec[i] <- 2
        else {
            power.n.max.i <- power.n.max[i]
            if (power.n.max.i < power.i) {
                n.vec[i] <- NA
                warning(paste("Error in search algorithm for element ", 
                  i, ".  Try increasing the argument 'n.max'", 
                  sep = ""))
            }
            else {
                delta.over.sigma.i <- delta.over.sigma[i]
                alpha.i <- alpha[i]
                n.vec[i] <- uniroot(fcn.for.root, lower = 2, 
                  upper = n.max, f.lower = power.i - power.2.i, 
                  f.upper = power.i - power.n.max.i, power = power.i, 
                  delta.over.sigma = delta.over.sigma.i, alpha = alpha.i, 
                  sample.type = sample.type, alternative = alternative, 
                  approx = approx, tol = tol, maxiter = maxiter)$root
            }
        }
    }
    if (n2.constrained) {
        n1 <- (n.vec * n2)/(2 * n2 - n.vec)
        if (any(index <- !is.finite(n1) | n1 < 0)) {
            n1[index] <- NA
            warning(paste("One or more constrained values of 'n2' is(are)", 
                "too small given the associated values of", "'power', 'delta.over.sigma', and 'alpha'"))
        }
        names(n1) <- names(n2) <- NULL
        if (round.up) 
            n1 <- ceiling(n1)
        ret.val <- list(n1 = n1, n2 = n2)
    }
    else {
        if (round.up) 
            n.vec <- ceiling(n.vec)
        ret.val <- n.vec
    }
    ret.val
}
