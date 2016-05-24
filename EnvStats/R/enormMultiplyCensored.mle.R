enormMultiplyCensored.mle <-
function (x, censored, N, cen.levels, K, c.vec, n.cen, censoring.side, 
    ci, ci.method = "profile.likelihood", ci.type, conf.level, 
    ci.sample.size = n.no.cen, pivot.statistic = "z") 
{
    if (censoring.side == "left") {
        x <- -x
        cen.levels <- -cen.levels
    }
    n.no.cen <- N - n.cen
    x.bar <- mean(x[!censored])
    s.mme <- sqrt((n.no.cen - 1)/n.no.cen) * sd(x[!censored])
    h.vec <- c.vec/n.no.cen
    fcn <- function(theta, x.bar, s.mme, h.vec, cen.levels, c.vec) {
        mu <- theta[1]
        sigma <- theta[2]
        zeta <- (cen.levels - mu)/sigma
        Q <- hnorm(zeta)
        con <- x.bar - mu
        eq1 <- con + sigma * sum(h.vec * Q)
        eq2 <- (s.mme^2 + con^2) - sigma^2 * (1 - sum(zeta * 
            h.vec * Q))
        eq1^2 + eq2^2
    }
    init.vec <- enormMultiplyCensored.qq.reg(x = x, censored = censored, 
        N = N, cen.levels = cen.levels, K = K, c.vec = c.vec, 
        n.cen = n.cen, censoring.side = "right", ci = FALSE)$parameters
    parameters <- nlminb(start = init.vec, objective = fcn, lower = c(-Inf, 
        .Machine$double.eps), x.bar = x.bar, s.mme = s.mme, h.vec = h.vec, 
        cen.levels = cen.levels, c.vec = c.vec)$par
    names(parameters) <- c("mean", "sd")
    ret.params <- parameters
    if (censoring.side == "left") 
        ret.params["mean"] <- -ret.params["mean"]
    ret.list <- list(parameters = ret.params)
    if (ci) {
        ci.method <- match.arg(ci.method, c("normal.approx", 
            "profile.likelihood"))
        pivot.statistic <- match.arg(pivot.statistic, c("z", 
            "t"))
        muhat <- parameters[1]
        sdhat <- parameters[2]
        s2 <- sdhat^2
        zeta <- (cen.levels - muhat)/sdhat
        Q <- hnorm(zeta)
        Qp <- Q * (Q - zeta)
        lambda <- zeta * Qp + Q
        eta <- zeta * (lambda + Q)
        FishInfMat <- matrix(0, 2, 2)
        con <- n.no.cen/s2
        FishInfMat[1, 1] <- con * (1 + sum(h.vec * Qp))
        FishInfMat[1, 2] <- con * ((2 * (x.bar - muhat))/sdhat + 
            sum(h.vec * lambda))
        FishInfMat[2, 1] <- FishInfMat[1, 2]
        FishInfMat[2, 2] <- con * ((3 * (s.mme^2 + (x.bar - muhat)^2))/s2 - 
            1 + sum(h.vec * eta))
        var.cov.params <- solve(FishInfMat)
        dimnames(var.cov.params) <- list(c("mean", "sd"), c("mean", 
            "sd"))
        ci.type.arg <- ci.type
        if (ci.type != "two-sided" && censoring.side == "left") 
            ci.type.arg <- ifelse(ci.type == "upper", "lower", 
                "upper")
        ci.obj <- ci.normal.approx(theta.hat = muhat, sd.theta.hat = sqrt(var.cov.params[1, 
            1]), n = ci.sample.size, df = ci.sample.size - 1, 
            ci.type = ci.type.arg, alpha = 1 - conf.level, test.statistic = pivot.statistic)
        ci.obj$parameter <- "mean"
        if (censoring.side == "left") {
            var.cov.params[1, 2] <- -var.cov.params[1, 2]
            var.cov.params[2, 1] <- var.cov.params[1, 2]
            limits <- ci.obj$limits
            names.limits <- names(limits)
            ci.obj$limits <- -rev(limits)
            names(ci.obj$limits) <- names.limits
        }
        if (ci.method == "profile.likelihood") {
            if (censoring.side == "left") {
                x <- -x
            }
            loglik.at.mle <- loglikCensored(theta = ret.params, 
                x = x, censored = censored, censoring.side = censoring.side, 
                distribution = "norm")
            fcn.pl <- function(CL, loglik.at.mle, sd.mle, x, 
                censored, censoring.side, conf.level) {
                sd.mle.at.CL <- enormCensored.sd.mle.at.fixed.mean(fixed.mean = CL, 
                  sd.mle = sd.mle, x = x, censored = censored, 
                  censoring.side = censoring.side)
                (2 * (loglik.at.mle - loglikCensored(theta = c(CL, 
                  sd.mle.at.CL), x = x, censored = censored, 
                  censoring.side = censoring.side, distribution = "norm")) - 
                  qchisq(conf.level, df = 1))^2
            }
            limits <- ci.obj$limits
            names(limits) <- NULL
            switch(ci.type, `two-sided` = {
                LCL <- nlminb(start = limits[1], objective = fcn.pl, 
                  upper = ret.params["mean"], loglik.at.mle = loglik.at.mle, 
                  sd.mle = ret.params["sd"], x = x, censored = censored, 
                  censoring.side = censoring.side, conf.level = conf.level)$par
                UCL <- nlminb(start = limits[2], objective = fcn.pl, 
                  lower = ret.params["mean"], loglik.at.mle = loglik.at.mle, 
                  sd.mle = ret.params["sd"], x = x, censored = censored, 
                  censoring.side = censoring.side, conf.level = conf.level)$par
            }, lower = {
                LCL <- nlminb(start = limits[1], objective = fcn.pl, 
                  upper = ret.params["mean"], loglik.at.mle = loglik.at.mle, 
                  sd.mle = ret.params["sd"], x = x, censored = censored, 
                  censoring.side = censoring.side, conf.level = 1 - 
                    2 * (1 - conf.level))$par
                UCL = Inf
            }, upper = {
                LCL = -Inf
                UCL <- nlminb(start = limits[2], objective = fcn.pl, 
                  lower = ret.params["mean"], loglik.at.mle = loglik.at.mle, 
                  sd.mle = ret.params["sd"], x = x, censored = censored, 
                  censoring.side = censoring.side, conf.level = 1 - 
                    2 * (1 - conf.level))$par
            })
            names(LCL) <- "LCL"
            names(UCL) <- "UCL"
            ci.obj <- list(name = "Confidence", parameter = "mean", 
                limits = c(LCL, UCL), type = ifelse(ci.type == 
                  "two.sided", "two-sided", ci.type), method = "Profile Likelihood", 
                conf.level = conf.level)
        }
        oldClass(ci.obj) <- "intervalEstimateCensored"
        ret.list <- c(ret.list, list(var.cov.params = var.cov.params, 
            ci.obj = ci.obj))
    }
    ret.list
}
