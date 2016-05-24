elnormAltMultiplyCensored.mle <-
function (x, censored, N, cen.levels, K, c.vec, n.cen, censoring.side, 
    ci, ci.method = c("profile.likelihood", "cox", "delta"), 
    ci.type, conf.level, ci.sample.size = N - n.cen, pivot.statistic = c("z", 
        "t")) 
{
    enorm.list <- enormMultiplyCensored.mle(x = log(x), censored = censored, 
        N = N, cen.levels = log(cen.levels), K = K, c.vec = c.vec, 
        n.cen = n.cen, censoring.side = censoring.side, ci = ci, 
        ci.method = "normal.approx", ci.type = ci.type, conf.level = conf.level)
    log.parameters <- enorm.list$parameters
    meanlog <- log.parameters[1]
    sdlog <- log.parameters[2]
    mean <- exp(meanlog + sdlog^2/2)
    cv <- sqrt(exp(sdlog^2) - 1)
    parameters <- c(mean, cv)
    names(parameters) <- c("mean", "cv")
    ret.list <- list(parameters = parameters)
    if (ci) {
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
        }
        if (ci.method == "profile.likelihood") {
            loglik.at.mle <- loglikCensored(theta = ret.list$parameters, 
                x = x, censored = censored, censoring.side = censoring.side, 
                distribution = "lnormAlt")
            fcn <- function(CL, loglik.at.mle, mean.mle, cv.mle, 
                x, censored, censoring.side, conf.level) {
                cv.mle.at.CL <- elnormAltCensored.cv.mle.at.fixed.mean(fixed.mean = CL, 
                  mean.mle = mean.mle, cv.mle = cv.mle, x = x, 
                  censored = censored, censoring.side = censoring.side)
                (2 * (loglik.at.mle - loglikCensored(theta = c(CL, 
                  cv.mle.at.CL), x = x, censored = censored, 
                  censoring.side = censoring.side, distribution = "lnormAlt")) - 
                  qchisq(conf.level, df = 1))^2
            }
            limits <- ci.obj$limits
            names(limits) <- NULL
            switch(ci.type, `two-sided` = {
                LCL <- nlminb(start = limits[1], objective = fcn, 
                  lower = .Machine$double.eps, upper = mean, 
                  loglik.at.mle = loglik.at.mle, mean.mle = mean, 
                  cv.mle = cv, x = x, censored = censored, censoring.side = censoring.side, 
                  conf.level = conf.level)$par
                UCL <- nlminb(start = limits[2], objective = fcn, 
                  lower = mean, loglik.at.mle = loglik.at.mle, 
                  mean.mle = mean, cv.mle = cv, x = x, censored = censored, 
                  censoring.side = censoring.side, conf.level = conf.level)$par
            }, lower = {
                LCL <- nlminb(start = limits[1], objective = fcn, 
                  lower = .Machine$double.eps, upper = mean, 
                  loglik.at.mle = loglik.at.mle, mean.mle = mean, 
                  cv.mle = cv, x = x, censored = censored, censoring.side = censoring.side, 
                  conf.level = 1 - 2 * (1 - conf.level))$par
                UCL <- Inf
            }, upper = {
                LCL <- 0
                UCL <- nlminb(start = limits[2], objective = fcn, 
                  lower = mean, loglik.at.mle = loglik.at.mle, 
                  mean.mle = mean, cv.mle = cv, x = x, censored = censored, 
                  censoring.side = censoring.side, conf.level = 1 - 
                    2 * (1 - conf.level))$par
            })
            names(LCL) <- "LCL"
            names(UCL) <- "UCL"
            ci.obj$limits <- c(LCL, UCL)
        }
        oldClass(ci.obj) <- "intervalEstimateCensored"
        ci.obj$method <- switch(ci.method, delta = "Delta", cox = "Cox", 
            profile.likelihood = "Profile Likelihood")
        ret.list <- c(ret.list, list(ci.obj = ci.obj))
    }
    ret.list
}
