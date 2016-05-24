enormMultiplyCensored.impute.w.qq.reg <-
function (x, censored, N, cen.levels, K, c.vec, n.cen, censoring.side, 
    ci, ci.method = "normal.approx", ci.type, conf.level, prob.method = c("hirsch-stedinger", 
        "michael-schucany", "kaplan-meier", "nelson"), plot.pos.con = 0.375, 
    ci.sample.size = N - n.cen, lb.impute = -Inf, ub.impute = Inf, 
    pivot.statistic = c("z", "t")) 
{
    prob.method <- match.arg(prob.method)
    parameters <- enormMultiplyCensored.qq.reg(x = x, censored = censored, 
        N = N, , cen.levels = cen.levels, K = K, c.vec = c.vec, 
        n.cen = n.cen, censoring.side = censoring.side, ci = FALSE, 
        prob.method = prob.method, plot.pos.con = plot.pos.con)$parameters
    ppoints.list <- ppointsCensored(x = x, censored = censored, 
        censoring.side = censoring.side, prob.method = prob.method, 
        plot.pos.con = plot.pos.con)
    x.obs <- x[!censored]
    E.norm <- qnorm(ppoints.list$Cumulative.Probabilities[ppoints.list$Censored])
    x.impute <- E.norm * parameters[2] + parameters[1]
    if (any(index <- x.impute < lb.impute)) 
        x.impute[index] <- lb.impute
    if (any(index <- x.impute > ub.impute)) 
        x.impute[index] <- ub.impute
    new.x <- c(x.impute, x.obs)
    parameters[1] <- mean(new.x)
    parameters[2] <- sd(new.x)
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
