elnormAltMultiplyCensored.half.cen.level <-
function (x, censored, N, cen.levels, K, c.vec, n.cen, censoring.side, 
    ci, ci.method = "normal.approx", ci.type, conf.level, ci.sample.size = N - 
        n.cen, pivot.statistic = c("z", "t")) 
{
    x.cen <- x[censored]
    for (j in 1:K) {
        Tj <- cen.levels[j]
        x.cen[x.cen == Tj] <- Tj/2
    }
    x[censored] <- x.cen
    muhat <- mean(x)
    sdhat <- sd(x)
    parameters <- c(mean = muhat, cv = sdhat/muhat)
    if (ci) {
        ci.method <- match.arg(ci.method)
        pivot.statistic <- match.arg(pivot.statistic)
        ci.obj <- ci.normal.approx(theta.hat = muhat, sd.theta.hat = sdhat/sqrt(ci.sample.size), 
            n = ci.sample.size, df = ci.sample.size - 1, ci.type = ci.type, 
            alpha = 1 - conf.level, lb = 0, test.statistic = pivot.statistic)
        ci.obj$parameter <- "mean"
        return(list(parameters = parameters, ci.obj = ci.obj))
    }
    else return(list(parameters = parameters))
}
