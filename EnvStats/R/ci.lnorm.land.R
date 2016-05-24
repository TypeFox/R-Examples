ci.lnorm.land <-
function (meanlog, sdlog, n, ci.type = c("two-sided", "lower", 
    "upper"), alpha = 0.05) 
{
    ci.type <- match.arg(ci.type)
    ret.obj <- ci.land(lambda = 1/2, mu.hat = meanlog, sig.sq.hat = sdlog^2, 
        n, nu = n - 1, gamma.sq = n, ci.type = ci.type, conf.level = 1 - 
            alpha)
    ret.obj$limits <- exp(ret.obj$limits)
    ret.obj$parameter <- "mean"
    ret.obj
}
