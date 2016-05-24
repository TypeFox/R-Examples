tTestAlpha <-
function (n.or.n1, n2 = n.or.n1, delta.over.sigma = 0, power = 0.95, 
    sample.type = ifelse(!missing(n2) && !is.null(n2), "two.sample", 
        "one.sample"), alternative = "two.sided", approx = FALSE, 
    tol = 1e-07, maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(delta.over.sigma, 
        mode = "numeric") || !is.vector(power, mode = "numeric")) 
        stop("'n.or.n1', 'delta.over.sigma', and 'power' must be numeric vectors.")
    if (!all(is.finite(n.or.n1)) || !all(is.finite(delta.over.sigma)) || 
        !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.or.n1', 'delta.over.sigma', or 'power'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2")
    if (any(power <= 0) || any(power >= 1)) 
        stop("All values of 'power' must be greater than 0 and less than 1")
    if (sample.type == "two.sample" && !missing(n2)) {
        if (is.null(n2) || !is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (!all(is.finite(n2))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in 'n2'"))
        if (any(n2 < 2)) 
            stop("All values of 'n2' must be greater than or equal to 2")
    }
    arg.mat <- cbind.no.warn(n.or.n1 = as.vector(n.or.n1), n2 = as.vector(n2), 
        power = as.vector(power), delta.over.sigma = as.vector(delta.over.sigma))
    for (i in c("n.or.n1", "n2", "power", "delta.over.sigma")) assign(i, 
        arg.mat[, i])
    N <- nrow(arg.mat)
    alpha.vec <- rep(as.numeric(NA), N)
    fcn.for.root <- function(alpha, n.or.n1, n2, delta.over.sigma, 
        power, sample.type, alternative, approx) {
        power - tTestPower(n.or.n1 = n.or.n1, n2 = n2, delta.over.sigma = delta.over.sigma, 
            alpha = alpha, sample.type = sample.type, alternative = alternative, 
            approx = approx)
    }
    lower <- .Machine$double.eps
    upper <- 1 - .Machine$double.eps
    for (i in 1:N) {
        n.or.n1.i <- n.or.n1[i]
        n2.i <- n2[i]
        power.i <- power[i]
        delta.over.sigma.i <- delta.over.sigma[i]
        f.lower.i <- fcn.for.root(lower, n.or.n1 = n.or.n1.i, 
            n2 = n2.i, delta.over.sigma = delta.over.sigma.i, 
            power = power.i, sample.type = sample.type, alternative = alternative, 
            approx = approx)
        f.upper.i <- fcn.for.root(upper, n.or.n1 = n.or.n1.i, 
            n2 = n2.i, delta.over.sigma = delta.over.sigma.i, 
            power = power.i, sample.type = sample.type, alternative = alternative, 
            approx = approx)
        alpha.vec[i] <- uniroot(fcn.for.root, n.or.n1 = n.or.n1.i, 
            n2 = n2.i, delta.over.sigma = delta.over.sigma.i, 
            power = power.i, sample.type = sample.type, alternative = alternative, 
            approx = approx, lower = lower, upper = upper, f.lower = f.lower.i, 
            f.upper = f.upper.i, tol = tol, maxiter = maxiter)$root
    }
    alpha.vec
}
