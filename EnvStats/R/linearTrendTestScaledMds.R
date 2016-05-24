linearTrendTestScaledMds <-
function (n, x = lapply(n, seq), alpha = 0.05, power = 0.95, 
    alternative = "two.sided", two.sided.direction = "greater", 
    approx = FALSE, tol = 1e-07, maxiter = 1000) 
{
    if (missing(n) && missing(x)) 
        stop("You must supply either 'n' or 'x'")
    if (!missing(x)) {
        if (is.vector(x) && !is.list(x)) 
            x <- list(x)
        if (!is.list(x) || !all(sapply(x, function(i) {
            is.vector(i, mode = "numeric")
        }))) 
            stop(paste("'x' must be either a numeric vector or", 
                "a list in which each component of 'x'", "is a numeric vector"))
        if (any(sapply(x, function(i) {
            !all(is.finite(i))
        }))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in", 
                "any components of 'x'"))
        n <- sapply(x, length)
        n.unique <- sapply(x, function(y) length(unique(y)))
        if (any(n < 3) || any(n.unique < 2)) 
            stop(paste("All components of 'x' must contain at least 3 elements", 
                "and at least 2 distinct values."))
    }
    else {
        if (!is.vector(n, mode = "numeric")) 
            stop("'n' must be a numeric vector")
        if (any(is.na(n))) 
            stop(paste("Missing (NA) and Undefined (Nan) values", 
                "are not allowed in 'n'"))
        if (any(n < 3)) 
            stop("All values of 'n' must be greater than or equal to 3.")
    }
    if (!is.vector(alpha, mode = "numeric") || !is.vector(power, 
        mode = "numeric")) 
        stop("'alpha', and 'power' must be numeric vectors.")
    if (!all(is.finite(alpha)) || !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'alpha', or 'power'"))
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power < alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than or equal to", 
            "the corresponding elements of 'alpha', and less than 1"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    alt.fac <- ifelse(alternative == "two.sided", 2, 1)
    df <- n - 2
    oossx <- 1/sqrt(sapply(x, function(x) (length(x) - 1) * var(x)))
    slope.over.sigma.vec <- oossx * (qt(1 - alpha/alt.fac, df) + 
        qt(power, df))
    index <- power == alpha
    slope.over.sigma.vec[index] <- 0
    if (!approx) {
        alt <- ifelse(alternative == "less", "greater", alternative)
        arg.mat <- cbind.no.warn(n = as.vector(n), list.index = 1:length(x), 
            power = as.vector(power), alpha = as.vector(alpha))
        n <- arg.mat[, "n"]
        list.index <- arg.mat[, "list.index"]
        power <- arg.mat[, "power"]
        alpha <- arg.mat[, "alpha"]
        x <- x[list.index]
        N <- nrow(arg.mat)
        fcn.for.root <- function(slope.over.sigma, x, power, 
            alpha, alternative, approx) {
            power - linearTrendTestPower(x = x, slope.over.sigma = slope.over.sigma, 
                alpha = alpha, alternative = alternative, approx = approx)
        }
        for (i in (1:N)[!index]) {
            x.i <- x[i]
            power.i <- power[i]
            alpha.i <- alpha[i]
            slope.over.sigma.i <- slope.over.sigma.vec[i]
            upper <- 2 * slope.over.sigma.i
            power.upper <- linearTrendTestPower(x = x.i, slope.over.sigma = upper, 
                alpha = alpha.i, alternative = alt, approx = FALSE)
            upper.too.small <- power.upper <= power.i
            iter <- 1
            while (upper.too.small && iter <= maxiter) {
                upper <- 2 * upper
                power.upper <- linearTrendTestPower(x = x.i, 
                  slope.over.sigma = upper, alpha = alpha.i, 
                  alternative = alt, approx = FALSE)
                upper.too.small <- power.upper <= power.i
                iter <- iter + 1
            }
            if (iter > maxiter) 
                stop("Error in search algorithm.  Try increasing the value of the argument 'maxiter'")
            slope.over.sigma.vec[i] <- uniroot(fcn.for.root, 
                lower = 0, upper = upper, f.lower = power.i - 
                  alpha.i, f.upper = power.i - power.upper, x = x.i, 
                power = power.i, alpha = alpha.i, alternative = alt, 
                approx = FALSE, tol = tol, maxiter = maxiter)$root
        }
    }
    if (alternative == "less" || (alternative == "two.sided" && 
        two.sided.direction == "less")) 
        slope.over.sigma.vec <- -slope.over.sigma.vec
    slope.over.sigma.vec
}
