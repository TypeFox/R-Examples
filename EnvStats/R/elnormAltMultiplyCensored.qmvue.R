elnormAltMultiplyCensored.qmvue <-
function (x, censored, N, cen.levels, K, c.vec, n.cen, censoring.side, 
    ci, ci.method = c("delta", "cox"), ci.type, conf.level, ci.sample.size = N - 
        n.cen, pivot.statistic = c("z", "t")) 
{
    enorm.list <- enormMultiplyCensored.mle(x = log(x), censored = censored, 
        N = N, cen.levels = log(cen.levels), K = K, c.vec = c.vec, 
        n.cen = n.cen, censoring.side = censoring.side, ci = ci, 
        ci.method = "normal.approx", ci.type = ci.type, conf.level = conf.level)
    log.parameters <- enorm.list$parameters
    meanlog <- log.parameters[1]
    sdlog <- log.parameters[2]
    s2 <- sdlog^2
    mean <- exp(meanlog) * finneys.g(N - 1, s2/2)
    sd <- sqrt(exp(2 * meanlog) * (finneys.g(N - 1, 2 * s2) - 
        finneys.g(N - 1, (s2 * (N - 2))/(N - 1))))
    cv <- sd/mean
    parameters <- c(mean, cv)
    names(parameters) <- c("mean", "cv")
    if (ci) {
        sep.string <- paste("\n", space(33), sep = "")
        ci.method <- match.arg(ci.method)
        pivot.statistic <- match.arg(pivot.statistic)
        V <- enorm.list$var.cov.params
        if (ci.method == "delta") {
            lambda.vec <- c(mean, sdlog * mean)
            var.mean <- lambda.vec %*% V %*% lambda.vec
            ci.obj <- ci.normal.approx(theta.hat = mean, sd.theta.hat = sqrt(var.mean), 
                n = ci.sample.size, df = ci.sample.size - 1, 
                ci.type = ci.type, alpha = 1 - conf.level, lb = 0, 
                test.statistic = pivot.statistic)
            ci.obj$parameter <- "mean"
            ci.obj$method <- paste(ci.obj$method, "Based on Delta Method", 
                sep = sep.string)
        }
        else {
            beta <- log(mean)
            sd.beta <- sqrt(V[1, 1] + 2 * sdlog * V[1, 2] + sdlog^2 * 
                V[2, 2])
            ci.obj <- ci.normal.approx(theta.hat = beta, sd.theta.hat = sd.beta, 
                n = ci.sample.size, df = ci.sample.size - 1, 
                ci.type = ci.type, alpha = 1 - conf.level, test.statistic = pivot.statistic)
            ci.obj$limits <- exp(ci.obj$limits)
            ci.obj$parameter <- "mean"
            ci.obj$method <- paste(ci.obj$method, "Based on Cox's Method", 
                sep = sep.string)
        }
        return(list(parameters = parameters, ci.obj = ci.obj))
    }
    else return(list(parameters = parameters))
}
