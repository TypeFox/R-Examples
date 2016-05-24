elnormAltSinglyCensored.bcmle <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = c("delta", 
    "cox"), ci.type, conf.level, ci.sample.size = En, pivot.statistic = c("z", 
    "t")) 
{
    pivot.statistic <- match.arg(pivot.statistic)
    enorm.list <- enormSinglyCensored.mle(x = log(x), censored = censored, 
        N = N, T1 = log(T1), n.cen = n.cen, censoring.side = censoring.side, 
        ci = TRUE, ci.method = "normal.approx", ci.type = ci.type, 
        conf.level = conf.level, pivot.statistic = pivot.statistic)
    log.parameters <- enorm.list$parameters
    meanlog <- log.parameters[1]
    sdlog <- log.parameters[2]
    V <- enorm.list$var.cov.params
    rel.bias <- exp(0.5 * (V[1, 1] + 2 * sdlog * V[1, 2] + sdlog^2 * 
        V[2, 2]))
    mean <- exp(meanlog + sdlog^2/2)/rel.bias
    cv <- sqrt(exp(sdlog^2) - 1)
    parameters <- c(mean, cv)
    names(parameters) <- c("mean", "cv")
    if (ci) {
        sep.string <- paste("\n", space(33), sep = "")
        ci.method <- match.arg(ci.method)
        En <- enorm.list$ci.obj$sample.size
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
            sd.beta <- sqrt(2 * log(rel.bias))
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
