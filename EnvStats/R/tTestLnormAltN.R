tTestLnormAltN <-
function (ratio.of.means, cv = 1, alpha = 0.05, power = 0.95, 
    sample.type = ifelse(!is.null(n2), "two.sample", "one.sample"), 
    alternative = "two.sided", approx = FALSE, n2 = NULL, round.up = TRUE, 
    n.max = 5000, tol = 1e-07, maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(ratio.of.means, mode = "numeric") || !is.vector(cv, 
        mode = "numeric") || !is.vector(alpha, mode = "numeric") || 
        !is.vector(power, mode = "numeric")) 
        stop("'ratio.of.means', 'cv', 'alpha', and 'power' must be numeric vectors.")
    if (!all(is.finite(ratio.of.means)) || !all(is.finite(cv)) || 
        !all(is.finite(alpha)) || !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'ratio.of.means', 'cv', 'alpha', or 'power'"))
    if (any(ratio.of.means <= 0)) 
        stop("All values of 'ratio.of.means' must be positive")
    if (any(cv <= 0)) 
        stop("All values of 'cv' must be positive")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power <= alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than or equal to", 
            "the corresponding elements of 'alpha' and less than 1"))
    if (alternative == "greater" && any(ratio.of.means <= 1)) 
        stop(paste("When alternative='greater',", "all values of 'ratio.of.means' must be greater than 1."))
    if (alternative == "less" && any(ratio.of.means >= 1 | ratio.of.means <= 
        0)) 
        stop(paste("When alternative='less',", "all values of 'ratio.of.means' must be greater than 0", 
            "and less than 1."))
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    delta.over.sigma <- log(ratio.of.means)/sqrt(log(cv^2 + 1))
    ret.val <- tTestN(delta.over.sigma = delta.over.sigma, alpha = alpha, 
        power = power, sample.type = sample.type, alternative = alternative, 
        approx = approx, n2 = n2, round.up = round.up, n.max = n.max, 
        tol = tol, maxiter = maxiter)
    ret.val
}
