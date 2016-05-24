elnormAltSinglyCensored.impute.w.qq.reg <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = "normal.approx", 
    ci.type, conf.level, plot.pos.con = 0.375, ci.sample.size = N - 
        n.cen, lb.impute = 0, ub.impute = Inf, pivot.statistic = c("z", 
        "t")) 
{
    log.parameters <- enormSinglyCensored.qq.reg(x = log(x), 
        censored = censored, N = N, T1 = log(T1), n.cen = n.cen, 
        censoring.side = censoring.side, ci = FALSE, plot.pos.con = plot.pos.con)$parameters
    x.obs <- x[!censored]
    index <- 1:n.cen
    if (censoring.side == "right") 
        index <- N - index + 1
    E.norm <- qnorm(ppoints(N, a = plot.pos.con))[index]
    x.impute <- exp(E.norm * log.parameters[2] + log.parameters[1])
    if (any(index <- x.impute < lb.impute)) 
        x.impute[index] <- lb.impute
    if (any(index <- x.impute > ub.impute)) 
        x.impute[index] <- ub.impute
    new.x <- c(x.impute, x.obs)
    muhat <- mean(new.x)
    sdhat <- sd(new.x)
    parameters <- c(mean = muhat, cv = sdhat/muhat)
    ret.list <- list(parameters = parameters, plot.pos.con = plot.pos.con)
    if (ci) {
        ci.method <- match.arg(ci.method)
        pivot.statistic <- match.arg(pivot.statistic)
        ci.obj <- ci.normal.approx(theta.hat = muhat, sd.theta.hat = sdhat/sqrt(ci.sample.size), 
            n = ci.sample.size, df = ci.sample.size - 1, ci.type = ci.type, 
            alpha = 1 - conf.level, lb = 0, test.statistic = pivot.statistic)
        ci.obj$parameter <- "mean"
        ret.list <- c(ret.list, list(ci.obj = ci.obj))
    }
    ret.list
}
