enormMultiplyCensored.qq.reg <-
function (x, censored, N, cen.levels, K, c.vec, n.cen, censoring.side, 
    ci, ci.method = "normal.approx", ci.type, conf.level, prob.method = c("hirsch-stedinger", 
        "michael-schucany", "kaplan-meier", "nelson"), plot.pos.con = 0.375, 
    ci.sample.size = N - n.cen, pivot.statistic = c("z", "t")) 
{
    prob.method <- match.arg(prob.method)
    qq.list <- qqPlotCensored(x = x, censored = censored, censoring.side = censoring.side, 
        prob.method = prob.method, plot.pos.con = plot.pos.con, 
        plot.it = FALSE)
    parameters <- lm(qq.list$y ~ qq.list$x, subset = !qq.list$Censored)$coef
    names(parameters) <- c("mean", "sd")
    ret.list <- list(parameters = parameters, prob.method = prob.method, 
        plot.pos.con = plot.pos.con)
    if (ci) {
        ci.method <- match.arg(ci.method)
        pivot.statistic <- match.arg(pivot.statistic)
        ci.obj <- ci.normal.approx(theta.hat = parameters[1], 
            sd.theta.hat = parameters[2]/sqrt(ci.sample.size), 
            n = ci.sample.size, df = ci.sample.size - 1, ci.type = ci.type, 
            alpha = 1 - conf.level, test.statistic = pivot.statistic)
        ci.obj$parameter <- "mean"
        ret.list <- c(ret.list, list(ci.obj = ci.obj))
    }
    ret.list
}
