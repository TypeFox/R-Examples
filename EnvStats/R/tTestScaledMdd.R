tTestScaledMdd <-
function (n.or.n1, n2 = n.or.n1, alpha = 0.05, power = 0.95, 
    sample.type = ifelse(!missing(n2) && !is.null(n2), "two.sample", 
        "one.sample"), alternative = "two.sided", two.sided.direction = "greater", 
    approx = FALSE, tol = 1e-07, maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(alpha, 
        mode = "numeric") || !is.vector(power, mode = "numeric")) 
        stop("'n.or.n1', 'alpha', and 'power' must be numeric vectors.")
    if (!all(is.finite(n.or.n1)) || !all(is.finite(alpha)) || 
        !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.or.n1', 'alpha', or 'power'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power < alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than or equal to", 
            "the corresponding elements of 'alpha', and less than 1"))
    if (sample.type == "two.sample" && !missing(n2)) {
        if (is.null(n2) || !is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (!all(is.finite(n2))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in 'n2'"))
        if (any(n2 < 2)) 
            stop("All values of 'n2' must be greater than or equal to 2")
    }
    alt.fac <- ifelse(alternative == "two.sided", 2, 1)
    if (sample.type == "two.sample") {
        df <- n.or.n1 + n2 - 2
        oorn <- 1/sqrt((n.or.n1 * n2)/(n.or.n1 + n2))
    }
    else {
        df <- n.or.n1 - 1
        oorn <- 1/sqrt(n.or.n1)
    }
    delta.over.sigma.vec <- oorn * (qt(1 - alpha/alt.fac, df) + 
        qt(power, df))
    index <- power == alpha
    delta.over.sigma.vec[index] <- 0
    if (!approx) {
        alt <- ifelse(alternative == "less", "greater", alternative)
        arg.mat <- cbind.no.warn(n.or.n1 = as.vector(n.or.n1), 
            n2 = as.vector(n2), power = as.vector(power), alpha = as.vector(alpha))
        for (i in c("n.or.n1", "n2", "power", "alpha")) assign(i, 
            arg.mat[, i])
        N <- nrow(arg.mat)
        fcn.for.root <- function(delta.over.sigma, n.or.n1, n2, 
            power, alpha, sample.type, alternative, approx) {
            power - tTestPower(n.or.n1 = n.or.n1, n2 = n2, delta.over.sigma = delta.over.sigma, 
                alpha = alpha, sample.type = sample.type, alternative = alternative, 
                approx = approx)
        }
        for (i in (1:N)[!index]) {
            n.or.n1.i <- n.or.n1[i]
            n2.i <- n2[i]
            power.i <- power[i]
            alpha.i <- alpha[i]
            delta.over.sigma.i <- delta.over.sigma.vec[i]
            upper <- 2 * delta.over.sigma.i
            power.upper <- tTestPower(n.or.n1 = n.or.n1.i, n2 = n2.i, 
                delta.over.sigma = upper, alpha = alpha.i, sample.type = sample.type, 
                alternative = alt, approx = FALSE)
            upper.too.small <- power.upper <= power.i
            iter <- 1
            while (upper.too.small && iter <= maxiter) {
                upper <- 2 * upper
                power.upper <- tTestPower(n.or.n1 = n.or.n1.i, 
                  n2 = n2.i, delta.over.sigma = upper, alpha = alpha.i, 
                  sample.type = sample.type, alternative = alt, 
                  approx = FALSE)
                upper.too.small <- power.upper <= power.i
                iter <- iter + 1
            }
            if (iter > maxiter) 
                stop("Error in search algorithm.  Try increasing the argument 'maxiter'")
            delta.over.sigma.vec[i] <- uniroot(fcn.for.root, 
                lower = 0, upper = upper, f.lower = power.i - 
                  alpha.i, f.upper = power.i - power.upper, n.or.n1 = n.or.n1.i, 
                n2 = n2.i, power = power.i, alpha = alpha.i, 
                sample.type = sample.type, alternative = alt, 
                approx = FALSE, tol = tol, maxiter = maxiter)$root
        }
    }
    if (alternative == "less" || (alternative == "two.sided" && 
        two.sided.direction == "less")) 
        delta.over.sigma.vec <- -delta.over.sigma.vec
    delta.over.sigma.vec
}
