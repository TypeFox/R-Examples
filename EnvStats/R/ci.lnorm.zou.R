ci.lnorm.zou <-
function (meanlog, sdlog, n, ci.type, alpha) 
{
    theta.2.hat <- sdlog^2/2
    pivot <- meanlog + theta.2.hat
    limits.mean <- ci.normal.approx(theta.hat = meanlog, sd.theta.hat = sdlog/sqrt(n), 
        n = n, df = n - 1, ci.type = ci.type, alpha = alpha, 
        test.statistic = "z")$limits
    limits.theta.2 <- ci.norm.var(sdhat = sdlog, n = n, ci.type = ci.type, 
        alpha = alpha)$limits/2
    switch(ci.type, `two-sided` = {
        lcl <- exp(pivot - sqrt((meanlog - limits.mean["LCL"])^2 + 
            (theta.2.hat - limits.theta.2["LCL"])^2))
        ucl <- exp(pivot + sqrt((limits.mean["UCL"] - meanlog)^2 + 
            (limits.theta.2["UCL"] - theta.2.hat)^2))
    }, lower = {
        lcl <- exp(pivot - sqrt((meanlog - limits.mean["LCL"])^2 + 
            (theta.2.hat - limits.theta.2["LCL"])^2))
        ucl <- Inf
    }, upper = {
        lcl <- -Inf
        ucl <- exp(pivot + sqrt((limits.mean["UCL"] - meanlog)^2 + 
            (limits.theta.2["UCL"] - theta.2.hat)^2))
    })
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = ci.type, method = "Zou", conf.level = 1 - 
            alpha, sample.size = n, dof = n - 1)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
