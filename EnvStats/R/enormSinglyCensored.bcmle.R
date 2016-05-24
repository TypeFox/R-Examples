enormSinglyCensored.bcmle <-
function (x, censored, N, T1, n.cen, censoring.side, ci, ci.method = c("normal.approx", 
    "normal.approx.w.cov"), ci.type, conf.level, ci.sample.size = En, 
    pivot.statistic = c("z", "t")) 
{
    parameters <- enormSinglyCensored.mle(x, censored, N, T1, 
        n.cen, censoring.side, ci = FALSE)$parameters
    muhat <- parameters["mean"]
    sdhat <- parameters["sd"]
    N0 <- N - n.cen
    Np1 <- N + 1
    p <- N0/Np1
    Bu <- -exp(2.692 - 5.439 * p)
    Bs <- -((0.312 + 0.859 * p)^(-2))
    muhat <- muhat - (ifelse(censoring.side == "left", -1, 1) * 
        sdhat * Bu)/Np1
    sdhat <- sdhat - (sdhat * Bs)/Np1
    parameters[] <- c(muhat, sdhat)
    if (ci) {
        ci.method <- match.arg(ci.method)
        zeta.hat <- (T1 - muhat)/sdhat
        eta.hat <- ifelse(censoring.side == "left", zeta.hat, 
            -zeta.hat)
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
        mu.12 <- con * (-phi.12/denom)
        if (censoring.side == "right") 
            mu.12 <- -mu.12
        mu.22 <- con * (phi.11/denom)
        con <- (sdhat^2)/En
        var.muhat <- con * mu.11
        var.sdhat <- con * mu.22
        cov.muhat.sdhat <- con * mu.12
        var.cov.params <- matrix(c(var.muhat, cov.muhat.sdhat, 
            cov.muhat.sdhat, var.sdhat), 2, 2)
        dimnames(var.cov.params) <- list(c("mean", "sd"), c("mean", 
            "sd"))
        if (ci.method == "normal.approx") {
            pivot.statistic <- match.arg(pivot.statistic)
            ci.obj <- ci.normal.approx(theta.hat = muhat, sd.theta.hat = sqrt(var.muhat), 
                n = ci.sample.size, df = ci.sample.size - 1, 
                ci.type = ci.type, alpha = 1 - conf.level, test.statistic = pivot.statistic)
            ci.obj$parameter <- "mean"
        }
        else ci.obj <- ci.norm.w.cov(muhat = muhat, sdhat = sdhat, 
            var.muhat = var.muhat, var.sdhat = var.sdhat, cov.muhat.sdhat = cov.muhat.sdhat, 
            n = ci.sample.size, df = ci.sample.size - 1, ci.type = ci.type, 
            alpha = 1 - conf.level)
        return(list(parameters = parameters, var.cov.params = var.cov.params, 
            ci.obj = ci.obj))
    }
    else return(list(parameters = parameters))
}
