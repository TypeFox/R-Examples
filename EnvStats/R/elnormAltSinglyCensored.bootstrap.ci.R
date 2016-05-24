elnormAltSinglyCensored.bootstrap.ci <-
function (x, censored, N, T1, censoring.side, est.fcn, ci.type, 
    conf.level, n.bootstraps, obs.mean, ...) 
{
    boot.vec <- numeric(n.bootstraps)
    too.few.obs.count <- 0
    no.cen.obs.count <- 0
    x.no.cen <- x[!censored]
    for (i in 1:n.bootstraps) {
        index <- sample(N, replace = TRUE)
        new.x <- x[index]
        new.censored <- censored[index]
        new.n.cen <- sum(new.censored)
        if ((N - new.n.cen) < 2) {
            too.few.obs.count <- too.few.obs.count + 1
            i <- i - 1
            next
        }
        if (new.n.cen == 0) {
            boot.vec[i] <- elnormAlt(new.x)$parameters[1]
            no.cen.obs.count <- no.cen.obs.count + 1
        }
        else boot.vec[i] <- do.call(est.fcn, list(x = new.x, 
            censored = new.censored, N = N, T1 = T1, n.cen = new.n.cen, 
            censoring.side = censoring.side, ci = FALSE, ...))$parameters[1]
    }
    alpha <- 1 - conf.level
    if (ci.type == "two.sided") 
        alpha <- alpha/2
    ci.limits.pct <- switch(ci.type, `two-sided` = quantile(boot.vec, 
        probs = c(alpha, 1 - alpha)), lower = c(quantile(boot.vec, 
        probs = alpha), Inf), upper = c(0, quantile(boot.vec, 
        probs = conf.level)))
    compute.bca <- length(unique(x.no.cen)) >= 3
    if (compute.bca) {
        za <- qnorm(alpha)
        z0 <- qnorm(sum(boot.vec <= obs.mean)/n.bootstraps)
        jack.vec <- elnormAltSinglyCensored.jackknife(x = x, 
            censored = censored, N = N, T1 = T1, censoring.side = censoring.side, 
            est.fcn = est.fcn, ci.type = ci.type, conf.level = conf.level, 
            ...)
        num <- sum(as.vector(scale(jack.vec, scale = FALSE))^3)
        denom <- 6 * (((length(jack.vec) - 1) * var(jack.vec))^(3/2))
        a <- num/denom
        ci.limits.bca <- switch(ci.type, `two-sided` = {
            alpha1 <- pnorm(z0 + (z0 + za)/(1 - a * (z0 + za)))
            alpha2 <- pnorm(z0 + (z0 - za)/(1 - a * (z0 - za)))
            quantile(boot.vec, probs = c(alpha1, alpha2))
        }, lower = {
            alpha1 <- pnorm(z0 + (z0 + za)/(1 - a * (z0 + za)))
            c(quantile(boot.vec, probs = alpha1), Inf)
        }, upper = {
            alpha2 <- pnorm(z0 + (z0 - za)/(1 - a * (z0 - za)))
            c(0, quantile(boot.vec, probs = alpha2))
        })
    }
    else ci.limits.bca <- switch(ci.type, `two-sided` = c(NA, 
        NA), lower = c(NA, Inf), upper = c(-Inf, NA))
    ci.limits <- c(ci.limits.pct, ci.limits.bca)
    names(ci.limits) <- c("Pct.LCL", "Pct.UCL", "BCa.LCL", "BCa.UCL")
    ret.obj <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = ci.type, method = "Bootstrap", 
        conf.level = conf.level, n.bootstraps = n.bootstraps, 
        too.few.obs.count = too.few.obs.count, no.cen.obs.count = no.cen.obs.count)
    oldClass(ret.obj) <- "intervalEstimateCensored"
    ret.obj
}
