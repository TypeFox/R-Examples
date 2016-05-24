ci.gamma.chisq.approx <-
function (x, shape, shape.est.method, ci.type, conf.level) 
{
    x.bar <- mean(x)
    n <- length(x)
    alpha <- 1 - conf.level
    switch(ci.type, `two-sided` = {
        alpha <- (1 - conf.level)/2
        LCL <- 2 * n * shape * x.bar/qchisq(p = 1 - alpha/2, 
            df = 2 * n * shape)
        UCL <- 2 * n * shape * x.bar/qchisq(p = alpha/2, df = 2 * 
            n * shape)
    }, lower = {
        LCL <- 2 * n * shape * x.bar/qchisq(p = 1 - alpha, df = 2 * 
            n * shape)
        UCL <- Inf
    }, upper = {
        LCL <- 0
        UCL <- 2 * n * shape * x.bar/qchisq(p = alpha, df = 2 * 
            n * shape)
    })
    ci.limits <- c(LCL, UCL)
    names(ci.limits) <- c("LCL", "UCL")
    interval <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = ci.type, method = paste("Chi-square approximation\n", 
            space(33), "using ", shape.est.method, " of 'shape'", 
            sep = ""), conf.level = conf.level, sample.size = length(x))
    oldClass(interval) <- "intervalEstimate"
    interval
}
