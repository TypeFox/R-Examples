egammaAltCensored.mle <-
function (x, censored, censoring.side, ci, ci.method = "profile.likelihood", 
    ci.type = "two-sided", conf.level, ci.sample.size = sum(!censored), 
    pivot.statistic = "z") 
{
    fcn <- function(theta, x, censored, censoring.side) {
        -loglikCensored(theta = theta, x = x, censored = censored, 
            censoring.side = censoring.side, distribution = "gammaAlt")
    }
    if (censoring.side == "left") {
        dum.x <- x
        dum.x[censored] <- dum.x[censored]/2
        par.init <- egammaAlt(x = dum.x, ci = FALSE)$parameters
    }
    else {
        par.init <- egammaAlt(x = x, ci = FALSE)$parameters
    }
    parameters <- nlminb(start = par.init, objective = fcn, x = x, 
        censored = censored, censoring.side = censoring.side, 
        lower = .Machine$double.eps)$par
    names(parameters) <- c("mean", "cv")
    if (ci) {
        opt.list <- optim(par = parameters, fn = fcn, x = x, 
            censored = censored, censoring.side = censoring.side, 
            hessian = ci)
        parameters <- opt.list$par
    }
    ret.list <- list(parameters = parameters)
    if (ci) {
        ci.method <- match.arg(ci.method, c("profile.likelihood", 
            "normal.approx"))
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        sd.mean.mle <- sqrt(solve(opt.list$hessian)["mean", "mean"])
        pivot.statistic <- match.arg(pivot.statistic, c("z", 
            "t"))
        ci.obj <- ci.normal.approx(theta.hat = parameters["mean"], 
            sd.theta.hat = sd.mean.mle, n = ci.sample.size, df = ci.sample.size - 
                1, ci.type = ci.type, alpha = 1 - conf.level, 
            test.statistic = pivot.statistic)
        ci.obj$parameter <- "mean"
        if (ci.method == "profile.likelihood") {
            limits <- ci.obj$limits
            names(limits) <- NULL
            ci.obj <- ci.gammaAltCensored.profile.likelihood(x = x, 
                censored = censored, censoring.side = censoring.side, 
                mean.mle = parameters["mean"], cv.mle = parameters["cv"], 
                ci.type = ci.type, conf.level = conf.level, LCL.normal.approx = limits[1], 
                UCL.normal.approx = limits[2])
        }
        ret.list <- c(ret.list, list(ci.obj = ci.obj))
    }
    ret.list
}
