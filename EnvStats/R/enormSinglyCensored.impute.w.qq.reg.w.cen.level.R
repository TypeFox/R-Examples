enormSinglyCensored.impute.w.qq.reg.w.cen.level <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = "normal.approx", 
    ci.type, conf.level, plot.pos.con = 0.375, ci.sample.size = N - 
        n.cen, pivot.statistic = c("z", "t"), lb.impute = -Inf, 
    ub.impute = Inf) 
{
    parameters <- enormSinglyCensored.qq.reg.w.cen.level(x, censored, 
        N, T1, n.cen, censoring.side, ci = FALSE, plot.pos.con = plot.pos.con)$parameters
    x.obs <- x[!censored]
    index <- 1:n.cen
    if (censoring.side == "right") 
        index <- N - index + 1
    E.norm <- qnorm(ppoints(N, a = plot.pos.con))[index]
    x.impute <- E.norm * parameters[2] + parameters[1]
    if (any(index <- x.impute < lb.impute)) 
        x.impute[index] <- lb.impute
    if (any(index <- x.impute > ub.impute)) 
        x.impute[index] <- ub.impute
    new.x <- c(x.impute, x.obs)
    parameters[1] <- mean(new.x)
    parameters[2] <- sd(new.x)
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
