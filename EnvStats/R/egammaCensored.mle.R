egammaCensored.mle <-
function (x, censored, censoring.side, ci, ci.method = "profile.likelihood", 
    ci.type, conf.level, ci.sample.size = sum(!censored), pivot.statistic = "z") 
{
    fcn <- function(theta, x, censored, censoring.side) {
        -loglikCensored(theta = theta, x = x, censored = censored, 
            censoring.side = censoring.side, distribution = "gamma")
    }
    if (censoring.side == "left") {
        dum.x <- x
        dum.x[censored] <- dum.x[censored]/2
        par.init <- egamma(x = dum.x, ci = FALSE)$parameters
    }
    else {
        par.init <- egamma(x = x, ci = FALSE)$parameters
    }
    parameters <- nlminb(start = par.init, objective = fcn, x = x, 
        censored = censored, censoring.side = censoring.side, 
        lower = 1e-07)$par
    names(parameters) <- c("shape", "scale")
    ret.list <- list(parameters = parameters)
    if (ci) {
        ci.method <- match.arg(ci.method, c("normal.approx", 
            "profile.likelihood"))
        pivot.statistic <- match.arg(pivot.statistic, c("z", 
            "t"))
        ci.obj <- egammaAltCensored.mle(x = x, censored = censored, 
            censoring.side = censoring.side, ci = ci, ci.method = ci.method, 
            ci.type = ci.type, conf.level = conf.level, ci.sample.size = ci.sample.size, 
            pivot.statistic = pivot.statistic)$ci.obj
        ret.list <- c(ret.list, list(ci.obj = ci.obj))
    }
    ret.list
}
