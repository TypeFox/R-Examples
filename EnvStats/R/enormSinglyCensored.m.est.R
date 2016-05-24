enormSinglyCensored.m.est <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = c("normal.approx", 
    "normal.approx.w.cov"), ci.type, conf.level, ci.sample.size = N - 
    n.cen, pivot.statistic = c("z", "t"), t.df = 3) 
{
    if (censoring.side == "right") {
        x <- -x
        T1 <- -T1
    }
    est.list <- M.singly.censored(x, T1, t.df = t.df)
    parameters <- est.list$parameters
    if (censoring.side == "right") {
        parameters[1] <- -parameters[1]
    }
    names(parameters) <- c("mean", "sd")
    if (ci) {
        muhat <- parameters[1]
        sdhat <- parameters[2]
        ci.method <- match.arg(ci.method)
        var.cov.params <- est.list$var.cov.params
        dimnames(var.cov.params) <- list(c("mean", "sd"), c("mean", 
            "sd"))
        var.muhat <- var.cov.params[1, 1]
        if (ci.method == "normal.approx") {
            pivot.statistic <- match.arg(pivot.statistic)
            ci.obj <- ci.normal.approx(theta.hat = muhat, sd.theta.hat = sqrt(var.muhat), 
                n = ci.sample.size, df = ci.sample.size - 1, 
                ci.type = ci.type, alpha = 1 - conf.level, test.statistic = pivot.statistic)
            ci.obj$parameter <- "mean"
        }
        else {
            var.sdhat <- var.cov.params[2, 2]
            cov.muhat.sdhat <- var.cov.params[1, 2]
            ci.obj <- ci.norm.w.cov(muhat = muhat, sdhat = sdhat, 
                var.muhat = var.muhat, var.sdhat = var.sdhat, 
                cov.muhat.sdhat = cov.muhat.sdhat, n = ci.sample.size, 
                df = ci.sample.size - 1, ci.type = ci.type, alpha = 1 - 
                  conf.level)
        }
        return(list(parameters = parameters, var.cov.params = var.cov.params, 
            ci.obj = ci.obj))
    }
    else return(list(parameters = parameters))
}
