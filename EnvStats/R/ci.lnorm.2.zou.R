ci.lnorm.2.zou <-
function (meanlog.x, sdlog.x, n.x, meanlog.y, sdlog.y, n.y, ci.type, 
    alpha) 
{
    Mhat.x <- exp(meanlog.x + sdlog.x^2/2)
    Mhat.y <- exp(meanlog.y + sdlog.y^2/2)
    delta.hat <- Mhat.x - Mhat.y
    limits.x <- ci.lnorm.zou(meanlog.x, sdlog.x, n.x, ci.type, 
        alpha)$limits
    limits.y <- ci.lnorm.zou(meanlog.y, sdlog.y, n.y, ci.type, 
        alpha)$limits
    lcl <- delta.hat - sqrt((Mhat.x - limits.x["LCL"])^2 + (limits.y["UCL"] - 
        Mhat.y)^2)
    ucl <- delta.hat + sqrt((limits.x["UCL"] - Mhat.x)^2 + (Mhat.y - 
        limits.y["LCL"])^2)
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = "mean.x - mean.y", 
        limits = ci.limits, type = ci.type, method = "Zou", conf.level = 1 - 
            alpha, sample.size = c(n.x = n.x, n.y = n.y), dof = c(dof.x = n.x - 
            1, dof.y = n.y - 1))
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
