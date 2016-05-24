ci.binom.adjusted.Wald <-
function (x, size, alpha = 0.05, ci.type = c("two-sided", "lower", 
    "upper")) 
{
    if (length(x) != 1 || length(size) != 1) 
        stop("x and size must be vectors of length 1")
    ci.type <- match.arg(ci.type)
    alpha.const <- ifelse(ci.type == "two-sided", alpha/2, alpha)
    const <- qnorm(1 - alpha.const)^2
    ci.obj <- ci.binom.Wald(x = x + const/2, size = size + const, 
        alpha = alpha, ci.type = ci.type, correct = FALSE, var.denom = "n", 
        warn = FALSE)
    method <- "Adjusted Wald normal approximation"
    ret.obj <- list(name = "Confidence", parameter = "prob", 
        limits = ci.obj$limits, type = ci.type, method = method, 
        conf.level = 1 - alpha, sample.size = size)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
