enormSinglyCensored.qq.reg.w.cen.level <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = "normal.approx", 
    ci.type, conf.level, plot.pos.con = 0.375, ci.sample.size = N - 
        n.cen, pivot.statistic = c("z", "t")) 
{
    qq.list <- qqPlot(x, plot.pos.con = plot.pos.con, plot.it = FALSE)
    if (censoring.side == "left") 
        index <- n.cen:N
    else index <- 1:(N - n.cen + 1)
    parameters <- lm(qq.list$y ~ qq.list$x, subset = index)$coef
    names(parameters) <- c("mean", "sd")
    ret.list <- list(parameters = parameters, plot.pos.con = plot.pos.con)
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
