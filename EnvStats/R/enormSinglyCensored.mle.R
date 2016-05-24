enormSinglyCensored.mle <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = "profile.likelihood", 
    ci.type = "two-sided", conf.level, ci.sample.size = En, pivot.statistic = "z") 
{
    if (censoring.side == "left") {
        x <- -x
        T1 <- -T1
    }
    n.no.cen <- N - n.cen
    x.bar <- mean(x[!censored])
    s.mme <- sqrt((n.no.cen - 1)/n.no.cen) * sd(x[!censored])
    fcn <- function(theta, x.bar, s.mme, N, T1, n.cen) {
        mu <- theta[1]
        sigma <- theta[2]
        zeta <- (T1 - mu)/sigma
        h <- n.cen/(N - n.cen)
        Q <- hnorm(zeta)
        con <- x.bar - mu
        eq1 <- con + sigma * h * Q
        eq2 <- (s.mme^2 + con^2) - sigma^2 * (1 - zeta * h * 
            Q)
        eq1^2 + eq2^2
    }
    init.vec <- enormSinglyCensored.qq.reg(x = x, censored = censored, 
        N = N, T1 = T1, n.cen = n.cen, censoring.side = "right", 
        ci = FALSE)$parameters
    parameters <- nlminb(start = init.vec, objective = fcn, lower = c(-Inf, 
        .Machine$double.eps), x.bar = x.bar, s.mme = s.mme, N = N, 
        T1 = T1, n.cen = n.cen)$par
    names(parameters) <- c("mean", "sd")
    ret.params <- parameters
    if (censoring.side == "left") 
        ret.params["mean"] <- -ret.params["mean"]
    ret.list <- list(parameters = ret.params)
    if (ci) {
        ci.method <- match.arg(ci.method, c("normal.approx", 
            "normal.approx.w.cov", "profile.likelihood"))
        muhat <- parameters[1]
        sdhat <- parameters[2]
        zeta.hat <- (T1 - muhat)/sdhat
        eta.hat <- -zeta.hat
        snorm.eta.hat <- snorm(eta.hat)
        En <- N * snorm.eta.hat
        Q1 <- hnorm(eta.hat)
        Q2 <- hnorm(-eta.hat)
        phi.11 <- 1 + Q1 * (Q2 + eta.hat)
        phi.12 <- Q1 * (1 + eta.hat * (Q2 + eta.hat))
        phi.22 <- 2 + eta.hat * phi.12
        denom <- phi.11 * phi.22 - (phi.12^2)
        con <- (1/snorm.eta.hat)
        mu.11 <- con * (phi.22/denom)
        mu.12 <- con * (phi.12/denom)
        mu.22 <- con * (phi.11/denom)
        con <- (sdhat^2)/En
        var.muhat <- con * mu.11
        var.sdhat <- con * mu.22
        cov.muhat.sdhat <- con * mu.12
        var.cov.params <- matrix(c(var.muhat, cov.muhat.sdhat, 
            cov.muhat.sdhat, var.sdhat), 2, 2)
        dimnames(var.cov.params) <- list(c("mean", "sd"), c("mean", 
            "sd"))
        ci.type.arg <- ci.type
        if (ci.type != "two-sided" && censoring.side == "left") 
            ci.type.arg <- ifelse(ci.type == "upper", "lower", 
                "upper")
        if (ci.method %in% c("normal.approx", "profile.likelihood")) {
            pivot.statistic <- match.arg(pivot.statistic, c("z", 
                "t"))
            ci.obj <- ci.normal.approx(theta.hat = muhat, sd.theta.hat = sqrt(var.muhat), 
                n = ci.sample.size, df = ci.sample.size - 1, 
                ci.type = ci.type.arg, alpha = 1 - conf.level, 
                test.statistic = pivot.statistic)
            ci.obj$parameter <- "mean"
        }
        else ci.obj <- ci.norm.w.cov(muhat = muhat, sdhat = sdhat, 
            var.muhat = var.muhat, var.sdhat = var.sdhat, cov.muhat.sdhat = cov.muhat.sdhat, 
            n = ci.sample.size, df = ci.sample.size - 1, ci.type = ci.type.arg, 
            alpha = 1 - conf.level)
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
                UCL <- Inf
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
