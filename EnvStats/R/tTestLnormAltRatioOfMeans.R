tTestLnormAltRatioOfMeans <-
function (n.or.n1, n2 = n.or.n1, cv = 1, alpha = 0.05, power = 0.95, 
    sample.type = ifelse(!missing(n2), "two.sample", "one.sample"), 
    alternative = "two.sided", two.sided.direction = "greater", 
    approx = FALSE, tol = 1e-07, maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(cv, 
        mode = "numeric") || !is.vector(alpha, mode = "numeric") || 
        !is.vector(power, mode = "numeric")) 
        stop("'n.or.n1', 'cv', 'alpha', and 'power' must be numeric vectors.")
    if (!all(is.finite(n.or.n1)) || !all(is.finite(cv)) || !all(is.finite(alpha)) || 
        !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.or.n1', 'cv', 'alpha', or 'power'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2")
    if (any(cv <= 0)) 
        stop("All values of 'cv' must be positive")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power < alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than or equal to", 
            "the corresponding elements of 'alpha', and less than 1"))
    if (sample.type == "two.sample" && !missing(n2)) {
        if (!is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (!all(is.finite(n2))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in 'n2'"))
        if (any(n2 < 2)) 
            stop("All values of 'n2' must be greater than or equal to 2")
    }
    delta.over.sigma <- tTestScaledMdd(n.or.n1 = n.or.n1, n2 = n2, 
        alpha = alpha, power = power, sample.type = sample.type, 
        alternative = alternative, two.sided.direction = two.sided.direction, 
        approx = approx, tol = tol, maxiter = maxiter)
    ratio.of.means <- exp(delta.over.sigma * sqrt(log(cv^2 + 
        1)))
    ratio.of.means
}
