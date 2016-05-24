enparCensored.bootstrap.ci <-
function (x, censored, censoring.side, correct.se, left.censored.min, 
    right.censored.max, est.fcn, ci.type, conf.level, n.bootstraps, 
    obs.mean, obs.se.mean) 
{
    N <- length(x)
    boot.vec.mean <- numeric(n.bootstraps)
    boot.vec.t <- numeric(n.bootstraps)
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
            mu.hat <- mean(new.x)
            boot.vec.mean[i] <- mu.hat
            boot.vec.t[i] <- sqrt(N) * (mu.hat - obs.mean)/sd(new.x)
            no.cen.obs.count <- no.cen.obs.count + 1
        }
        else {
            params <- do.call(est.fcn, list(x = new.x, censored = new.censored, 
                censoring.side = censoring.side, correct.se = correct.se, 
                left.censored.min = left.censored.min, right.censored.max = right.censored.max, 
                ci = FALSE))$parameters
            mu.hat <- params["mean"]
            boot.vec.mean[i] <- mu.hat
            boot.vec.t[i] <- (mu.hat - obs.mean)/params["se.mean"]
        }
    }
    alpha <- 1 - conf.level
    if (ci.type == "two.sided") 
        alpha <- alpha/2
    ci.limits.pct <- switch(ci.type, `two-sided` = quantile(boot.vec.mean, 
        probs = c(alpha, 1 - alpha)), lower = c(quantile(boot.vec.mean, 
        probs = alpha), Inf), upper = c(0, quantile(boot.vec.mean, 
        probs = conf.level)))
    compute.bca <- length(unique(x.no.cen)) >= 3
    if (compute.bca) {
        za <- qnorm(alpha)
        z0 <- qnorm(sum(boot.vec.mean <= obs.mean)/n.bootstraps)
        jack.vec <- enparCensored.jackknife(x = x, censored = censored, 
            censoring.side = censoring.side, correct.se = correct.se, 
            left.censored.min = left.censored.min, right.censored.max = right.censored.max, 
            est.fcn = est.fcn)
        num <- sum(as.vector(scale(jack.vec, scale = FALSE))^3)
        denom <- 6 * (((length(jack.vec) - 1) * var(jack.vec))^(3/2))
        a <- num/denom
        ci.limits.bca <- switch(ci.type, `two-sided` = {
            alpha1 <- pnorm(z0 + (z0 + za)/(1 - a * (z0 + za)))
            alpha2 <- pnorm(z0 + (z0 - za)/(1 - a * (z0 - za)))
            quantile(boot.vec.mean, probs = c(alpha1, alpha2))
        }, lower = {
            alpha1 <- pnorm(z0 + (z0 + za)/(1 - a * (z0 + za)))
            c(quantile(boot.vec.mean, probs = alpha1), Inf)
        }, upper = {
            alpha2 <- pnorm(z0 + (z0 - za)/(1 - a * (z0 - za)))
            c(0, quantile(boot.vec.mean, probs = alpha2))
        })
    }
    else ci.limits.bca <- switch(ci.type, `two-sided` = c(NA, 
        NA), lower = c(NA, Inf), upper = c(0, NA))
    ci.limits.t <- switch(ci.type, `two-sided` = {
        t.quantiles <- quantile(boot.vec.t, probs = c(1 - alpha, 
            alpha))
        c(obs.mean - t.quantiles[1] * obs.se.mean, obs.mean - 
            t.quantiles[2] * obs.se.mean)
    }, lower = {
        t.quantiles <- quantile(boot.vec.t, probs = 1 - alpha)
        c(obs.mean - t.quantiles * obs.se.mean, Inf)
    }, upper = {
        t.quantiles <- quantile(boot.vec.t, probs = alpha)
        c(0, obs.mean - t.quantiles * obs.se.mean)
    })
    ci.limits <- c(ci.limits.pct, ci.limits.bca, ci.limits.t)
    names(ci.limits) <- c("Pct.LCL", "Pct.UCL", "BCa.LCL", "BCa.UCL", 
        "t.LCL", "t.UCL")
    ret.obj <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = ci.type, method = "Bootstrap", 
        conf.level = conf.level, sample.size = N, n.bootstraps = n.bootstraps, 
        too.few.obs.count = too.few.obs.count, no.cen.obs.count = no.cen.obs.count)
    oldClass(ret.obj) <- "intervalEstimateCensored"
    ret.obj
}
