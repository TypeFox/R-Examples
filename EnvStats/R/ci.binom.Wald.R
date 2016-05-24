ci.binom.Wald <-
function (x, size, alpha = 0.05, ci.type = c("two-sided", "lower", 
    "upper"), correct = TRUE, var.denom = c("n", "n-1"), warn = TRUE) 
{
    if (length(x) != 1 || length(size) != 1) 
        stop("x and size must be vectors of length 1")
    ci.type <- match.arg(ci.type)
    var.denom <- match.arg(var.denom)
    p.hat <- x/size
    q.hat <- 1 - p.hat
    if (warn && any(size * c(p.hat, q.hat) < 5) || p.hat < 0.2 || 
        p.hat > 0.8) 
        warning("Wald normal approximation for confidence interval may not be very accurate")
    s <- sqrt((p.hat * q.hat)/ifelse(var.denom == "n-1", (size - 
        1), size))
    switch(ci.type, `two-sided` = {
        ci.hw <- (qnorm(1 - (alpha/2)) * s) + correct * (1/(2 * 
            size))
        lcl <- pmax(0, p.hat - ci.hw)
        ucl <- pmin(1, p.hat + ci.hw)
    }, lower = {
        ci.hw <- (qnorm(1 - alpha) * s) + correct * (1/(2 * size))
        lcl <- pmax(0, p.hat - ci.hw)
        ucl <- 1
    }, upper = {
        ci.hw <- (qnorm(1 - alpha) * s) + correct * (1/(2 * size))
        lcl <- 0
        ucl <- pmin(1, p.hat + ci.hw)
    })
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    method <- paste("Wald normal approximation using '", var.denom, 
        "'\n", space(33), "for variance denominator", sep = "")
    if (correct) 
        method <- paste(method, "\n", space(33), "(With continuity correction)", 
            sep = "")
    ret.obj <- list(name = "Confidence", parameter = "prob", 
        limits = ci.limits, type = ci.type, method = method, 
        conf.level = 1 - alpha, sample.size = size)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
