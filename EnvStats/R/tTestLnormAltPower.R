tTestLnormAltPower <-
function (n.or.n1, n2 = n.or.n1, ratio.of.means = 1, cv = 1, 
    alpha = 0.05, sample.type = ifelse(!missing(n2), "two.sample", 
        "one.sample"), alternative = "two.sided", approx = FALSE) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(ratio.of.means, 
        mode = "numeric") || !is.vector(cv, mode = "numeric") || 
        !is.vector(alpha, mode = "numeric")) 
        stop("'n.or.n1', 'ratio.of.means', 'cv', and 'alpha' must be numeric vectors.")
    if (any(is.na(c(n.or.n1, ratio.of.means, cv)))) 
        stop(paste("Missing (NA) and Undefined (Nan) values", 
            "are not allowed in 'n.or.n1', 'ratio.of.means', or 'cv'"))
    if (!all(is.finite(alpha))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'alpha'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2.")
    if (any(ratio.of.means <= 0)) 
        stop("All values of 'ratio.of.means' must be positive")
    if (any(cv <= 0)) 
        stop("All values of 'cv' must be positive")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be between 0 and 1.")
    if (sample.type == "two.sample" && !missing(n2)) {
        if (!is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (any(is.na(n2))) 
            stop(paste("Missing (NA) and Undefined (Nan) values", 
                "are not allowed in 'n2'"))
        if (any(n2 < 2)) 
            stop("All values of 'n2' must be greater than or equal to 2")
    }
    delta.over.sigma <- log(ratio.of.means)/sqrt(log(cv^2 + 1))
    if (sample.type == "one.sample") {
        power <- tTestPower(n.or.n1 = n.or.n1, delta.over.sigma = delta.over.sigma, 
            alpha = alpha, sample.type = "one.sample", alternative = alternative, 
            approx = approx)
    }
    else {
        power <- tTestPower(n.or.n1 = n.or.n1, n2 = n2, delta.over.sigma = delta.over.sigma, 
            alpha = alpha, sample.type = "two.sample", alternative = alternative, 
            approx = approx)
    }
    power
}
