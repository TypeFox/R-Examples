ci.lnorm3.likelihood.profile <-
function (threshold, meanlog, sdlog, x, x1, n, eta.lb, eta.ub, 
    ci.type, alpha, ci.parameter) 
{
    conf.level <- 1 - alpha
    df <- n - 2
    eta.lmle <- -log(x1 - threshold)
    fcn.for.root <- function(eta, eta.lmle, x.weird, x1, n.weird, 
        crit.val) {
        threshold <- x1 - exp(-eta)
        y <- log(x.weird - threshold)
        ny <- length(y)
        meanlog <- mean(y)
        sdlog <- sqrt((ny - 1)/ny) * sd(y)
        threshold.lmle <- x1 - exp(-eta.lmle)
        y.lmle <- log(x.weird - threshold.lmle)
        meanlog.lmle <- mean(y.lmle)
        sdlog.lmle <- sqrt((ny - 1)/ny) * sd(y.lmle)
        n.weird * (2 * (meanlog - meanlog.lmle) + 2 * log(sdlog/sdlog.lmle)) - 
            crit.val
    }
    switch(ci.type, `two-sided` = {
        crit.val <- qchisq(conf.level, 1)
        f.lower.lcl <- fcn.for.root(eta.lb, eta.lmle, x, x1, 
            n, crit.val)
        f.upper.lcl <- fcn.for.root(eta.lmle - .Machine$double.eps, 
            eta.lmle, x, x1, n, crit.val)
        f.lower.ucl <- fcn.for.root(eta.lmle + .Machine$double.eps, 
            eta.lmle, x, x1, n, crit.val)
        f.upper.ucl <- fcn.for.root(eta.ub, eta.lmle, x, x1, 
            n, crit.val)
        if (sign(f.lower.lcl) == sign(f.upper.lcl) || sign(f.lower.ucl) == 
            sign(f.upper.ucl)) {
            warning(paste("Unable to solve for confidence limits for 'threshold'", 
                "in the interval", "[ mean(x) - 'threshold.lb.sd' * sd(x), min(x) ). ", 
                "Try changing the value of 'threshold.lb.sd'"))
            lcl <- ucl <- NA
        } else {
            eta.lcl <- uniroot(fcn.for.root, lower = eta.lb, 
                upper = eta.lmle - .Machine$double.eps, f.lower = f.lower.lcl, 
                f.upper = f.upper.lcl, eta.lmle = eta.lmle, x.weird = x, 
                x1 = x1, n.weird = n, crit.val = crit.val)$root
            eta.ucl <- uniroot(fcn.for.root, lower = eta.lmle + 
                .Machine$double.eps, upper = eta.ub, f.lower = f.lower.ucl, 
                f.upper = f.upper.ucl, eta.lmle = eta.lmle, x.weird = x, 
                x1 = x1, n.weird = n, crit.val = crit.val)$root
            if (eta.lcl == eta.lb || eta.ucl == eta.ub) {
                warning(paste("Unable to solve for confidence limits for 'threshold'", 
                  "in the interval", "[ mean(x) - 'threshold.lb.sd' * sd(x), min(x) ). ", 
                  "Try changing the value of 'threshold.lb.sd'"))
                lcl <- ucl <- NA
            } else {
                lcl <- x1 - exp(-eta.lcl)
                ucl <- x1 - exp(-eta.ucl)
            }
        }
    }, lower = {
        crit.val <- qchisq(1 - 2 * alpha, 1)
        ucl <- x1
        f.lower.lcl <- fcn.for.root(eta.lb, eta.lmle, x, x1, 
            n, crit.val)
        f.upper.lcl <- fcn.for.root(eta.lmle - .Machine$double.eps, 
            eta.lmle, x, x1, n, crit.val)
        if (sign(f.lower.lcl) == sign(f.upper.lcl)) {
            warning(paste("Unable to find lower confidence limit for 'threshold'", 
                "in the interval", "[ mean(x) - 'threshold.lb.sd' * sd(x), 'threshold.hat' ). ", 
                "Try changing the value of 'threshold.lb.sd'"))
            lcl <- NA
        } else {
            eta.lcl <- uniroot(fcn.for.root, lower = eta.lb, 
                upper = eta.lmle - .Machine$double.eps, f.lower = f.lower.lcl, 
                f.upper = f.upper.lcl, eta.lmle = eta.lmle, x.weird = x, 
                x1 = x1, n.weird = n, crit.val = crit.val)$root
            if (eta.lcl == eta.lb) {
                warning(paste("Unable to find lower confidence limit for 'threshold'", 
                  "in the interval", "[ mean(x) - 'threshold.lb.sd' * sd(x), 'threshold.hat' ). ", 
                  "Try changing the value of 'threshold.lb.sd'"))
                lcl <- NA
            } else {
                lcl <- x1 - exp(-eta.lcl)
            }
        }
    }, upper = {
        crit.val <- qchisq(1 - 2 * alpha, 1)
        lcl <- -Inf
        f.lower.ucl <- fcn.for.root(eta.lmle + .Machine$double.eps, 
            eta.lmle, x, x1, n, crit.val)
        f.upper.ucl <- fcn.for.root(eta.ub, eta.lmle, x, x1, 
            n, crit.val)
        if (sign(f.lower.ucl) == sign(f.upper.ucl)) {
            warning(paste("Unable to find upper confidence limit for 'threshold'", 
                "in the interval", "( 'threshold.hat', min(x) )"))
            ucl <- NA
        } else {
            eta.ucl <- uniroot(fcn.for.root, lower = eta.lmle + 
                .Machine$double.eps, upper = eta.ub, f.lower = f.lower.ucl, 
                f.upper = f.upper.ucl, eta.lmle = eta.lmle, x.weird = x, 
                x1 = x1, n.weird = n, crit.val = crit.val)$root
            if (eta.ucl == eta.ub) {
                warning(paste("Unable to find upper confidence limit for 'threshold'", 
                  "in the interval", "( 'threshold.hat', min(x) )"))
                ucl <- NA
            } else {
                ucl <- x1 - exp(-eta.ucl)
            }
        }
    })
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    if (ci.parameter == "median") {
        ci.obj.meanlog <- ci.normal.approx(meanlog, sdlog/sqrt(n), 
            n = n, df = df, ci.type = ci.type, alpha = alpha)
        ci.limits <- ci.limits + exp(ci.obj.meanlog$limits)
    }
    ret.obj <- list(name = "Confidence", parameter = ci.parameter, 
        limits = ci.limits, type = ci.type, method = "Likelihood Profile", 
        conf.level = conf.level, sample.size = n, dof = df)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
