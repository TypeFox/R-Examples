enormSinglyCensored.iterative.impute.w.qq.reg <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = "normal.approx", 
    ci.type, conf.level, plot.pos.con = 0.375, ci.sample.size = N - 
        n.cen, pivot.statistic = c("z", "t"), lb.impute = -Inf, 
    ub.impute = Inf, tol = 1e-06, convergence = c("relative", 
        "absolute"), max.iter = 100) 
{
    parameters <- enormSinglyCensored.impute.w.qq.reg(x, censored, 
        N, T1, n.cen, censoring.side, ci = FALSE, plot.pos.con = plot.pos.con, 
        lb.impute = lb.impute, ub.impute = ub.impute)$parameters
    parameters.old <- parameters
    x.obs <- x[!censored]
    index <- 1:n.cen
    if (censoring.side == "right") 
        index <- N - index + 1
    E.norm <- evNormOrdStats(N, approximate = TRUE)[index]
    convergence <- match.arg(convergence)
    iter <- 1
    diff <- tol + 1
    while (max(abs(diff)) >= tol && iter <= max.iter) {
        x.impute <- E.norm * parameters[2] + parameters[1]
        if (any(index <- x.impute < lb.impute)) 
            x.impute[index] <- lb.impute
        if (any(index <- x.impute > ub.impute)) 
            x.impute[index] <- ub.impute
        new.x <- c(x.impute, x.obs)
        parameters[] <- c(mean(new.x), sd(new.x))
        diff <- (parameters - parameters.old)
        if (convergence == "relative") 
            diff <- diff/parameters.old
        parameters.old <- parameters
        iter <- iter + 1
    }
    if (iter > max.iter) 
        warning(paste("No comvergence after", max.iter, "iterations"))
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
