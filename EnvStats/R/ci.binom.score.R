ci.binom.score <-
function (x, size, alpha = 0.05, ci.type = c("two-sided", "lower", 
    "upper"), correct = TRUE) 
{
    if (length(x) != 1 || length(size) != 1) 
        stop("x and size must be vectors of length 1")
    ci.type <- match.arg(ci.type)
    alternative <- switch(ci.type, `two-sided` = "two.sided", 
        lower = "greater", upper = "less")
    o.warn <- options(warn = -1)
    ci.limits <- prop.test(x = x, n = size, alternative = alternative, 
        conf.level = 1 - alpha, correct = correct)$conf.int
    options(o.warn)
    attributes(ci.limits) <- NULL
    names(ci.limits) <- c("LCL", "UCL")
    method <- "Score normal approximation"
    if (correct) 
        method <- paste(method, "\n", space(33), "(With continuity correction)", 
            sep = "")
    ret.obj <- list(name = "Confidence", parameter = "prob", 
        limits = ci.limits, type = ci.type, method = method, 
        conf.level = 1 - alpha, sample.size = size)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
