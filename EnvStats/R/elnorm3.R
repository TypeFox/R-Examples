elnorm3 <-
function (x, method = "lmle", ci = FALSE, ci.parameter = "threshold", 
    ci.method = "avar", ci.type = "two-sided", conf.level = 0.95, 
    threshold.lb.sd = 100) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector.")
    data.name <- deparse(substitute(x))
    method <- match.arg(method, c("lmle", "mme", "mmue", "mmme", 
        "royston.skew", "zero.skew"))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    if (ci) {
        ci.parameter <- match.arg(ci.parameter, c("threshold", 
            "median"))
        ci.method <- match.arg(ci.method, c("avar", "likelihood.profile", 
            "skewness"))
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        if (ci.method == "likelihood.profile" && method != "lmle") 
            stop(paste("You must set method='lmle' to use", "ci.method='likelihood.profile'"))
        if (ci.method == "skewness" && method != "zero.skew") 
            stop(paste("You must set method='zero.skew' to use", 
                "ci.method='skewness'"))
    }
    if (length(threshold.lb.sd) != 1 || !is.vector(threshold.lb.sd, 
        mode = "numeric") || threshold.lb.sd <= 0) 
        stop("'threshold.lb.sd' must be a positive scalar")
    n <- length(x)
    if (n < 3 || length(unique(x)) < 2) 
        stop("'x' must contain at least 3 non-missing distinct values")
    mean.x <- mean(x)
    var.x <- var(x)
    sd.x <- sqrt(var.x)
    x1 <- min(x)
    cf <- (n - 1)/n
    if (method == "royston.skew") {
        xn <- max(x)
        med.x <- median(x)
        q <- (xn - med.x)/(med.x - x1)
        if (q <= 1) 
            stop("Royston's skewness index indicates left-skewed data")
        threshold <- (x1 * xn - med.x^2)/(x1 + xn - 2 * med.x)
        y <- log(x - threshold)
        meanlog <- mean(y)
        sdlog <- sd(y)
    }
    else {
        m2 <- cf * var.x
        sqrt.b1 <- skewness(x, method = "moment")
        if (sqrt.b1 <= 0) 
            stop(paste("The sample skew is not positive. ", "Admissible moment estimates do not exist"))
        b1 <- sqrt.b1^2
        t1 <- 1 + b1/2
        t2 <- sqrt(t1^2 - 1)
        omega <- (t1 + t2)^(1/3) + (t1 - t2)^(1/3) - 1
        varlog <- log(omega)
        sdlog <- sqrt(varlog)
        meanlog <- 0.5 * log(m2/(omega * (omega - 1)))
        threshold <- mean.x - exp(meanlog + varlog/2)
        if (method != "mme") {
            meanlog <- 0.5 * log(var.x/(omega * (omega - 1)))
            threshold <- mean.x - exp(meanlog + varlog/2)
        }
        if (method != "mme" && method != "mmue") {
            fcn.to.min <- function(omega, s2, x.bar, x1, EZ1.n) {
                (s2/((x.bar - x1)^2) - (omega * (omega - 1))/((sqrt(omega) - 
                  exp(sqrt(log(omega)) * EZ1.n))^2))^2
            }
            nlminb.list <- nlminb(start = omega, objective = fcn.to.min, 
                lower = 1 + sqrt(.Machine$double.eps), s2 = var.x, 
                x.bar = mean.x, x1 = x1, EZ1.n = evNormOrdStatsScalar(r = 1, 
                  n = n))
            if (nlminb.list$convergence != 0) {
                warning("Unable to solve for 'omega'")
                threshold <- meanlog <- sdlog <- NA
                ci <- FALSE
            }
            else {
                omega <- nlminb.list$par
                varlog <- log(omega)
                sdlog <- sqrt(varlog)
                meanlog <- 0.5 * log(var.x/(omega * (omega - 
                  1)))
                threshold <- mean.x - exp(meanlog + varlog/2)
            }
        }
    }
    switch(method, mme = {
    }, mmue = {
    }, mmme = {
    }, lmle = {
        threshold.lb <- mean.x - threshold.lb.sd * sd.x
        threshold.ub <- x1 - sqrt(.Machine$double.eps)
        eta.lb <- -log(x1 - threshold.lb)
        eta.ub <- -log(x1 - threshold.ub)
        eta <- -log(x1 - threshold)
        neg.2.ll <- function(eta, x.weird, x1, n.weird) {
            threshold <- x1 - exp(-eta)
            y <- log(x.weird - threshold)
            ny <- length(y)
            meanlog <- mean(y)
            sdlog <- sqrt((ny - 1)/ny) * sd(y)
            n.weird * (1 + log(2 * pi) + 2 * meanlog + 2 * log(sdlog))
        }
        nlminb.list <- nlminb(start = eta, objective = neg.2.ll, 
            lower = eta.lb, upper = eta.ub, x.weird = x, x1 = x1, 
            n.weird = n)
        threshold <- x1 - exp(-nlminb.list$par)
        if (threshold == threshold.lb || threshold == threshold.ub) {
            warning(paste("Unable to solve for 'threshold' in the interval", 
                "[ mean(x) -", threshold.lb.sd, "* sd(x), min(x) ). ", 
                "Try changing the value of 'threshold.lb.sd'"))
            threshold <- meanlog <- sdlog <- NA
            ci <- FALSE
        } else {
            fcn.for.root.lmle <- function(threshold, x.weird, 
                n.weird) {
                y <- x.weird - threshold
                ly <- log(y)
                sly <- sum(ly)
                sum(1/y) * (sly - sum(ly^2) + (sly^2)/n.weird) - 
                  n.weird * sum(ly/y)
            }
            con <- 0.25 * abs(threshold)
            threshold <- uniroot(fcn.for.root.lmle, lower = max(threshold - 
                con, threshold.lb), upper = min(threshold + con, 
                mean(c(threshold, threshold.ub))), tol = .Machine$double.eps, 
                x.weird = x, n.weird = n)$root
            y <- log(x - threshold)
            ny <- length(y)
            meanlog <- mean(y)
            sdlog <- sqrt((ny - 1)/ny) * sd(y)
        }
    }, zero.skew = {
        threshold.lb <- mean.x - threshold.lb.sd * sd.x
        threshold.ub <- x1 - sqrt(.Machine$double.eps)
        fcn.for.root.zs <- function(threshold, x.weird) {
            y <- log(x.weird - threshold)
            sum((y - mean(y))^3)
        }
        uniroot.list <- uniroot(fcn.for.root.zs, lower = threshold.lb, 
            upper = threshold.ub, tol = .Machine$double.eps, 
            x.weird = x)
        threshold <- uniroot.list$root
        if (uniroot.list$f.root > sqrt(.Machine$double.eps) || 
            threshold == threshold.lb || threshold == threshold.ub) {
            warning(paste("Unable to solve for 'threshold' in the interval", 
                "[ mean(x) -", threshold.lb.sd, "* sd(x), min(x) ). ", 
                "Try changing the value of threshold.lb.sd"))
            threshold <- meanlog <- sdlog <- NA
            ci <- FALSE
        } else {
            y <- log(x - threshold)
            meanlog <- mean(y)
            sdlog <- sd(y)
        }
    }, royston.skew = {
    })
    ret.list <- list(distribution = "3-Parameter Lognormal", 
        sample.size = n, parameters = c(meanlog = meanlog, sdlog = sdlog, 
            threshold = threshold), n.param.est = 3, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        alpha <- 1 - conf.level
        beta <- exp(meanlog)
        beta2 <- beta^2
        switch(ci.method, avar = {
            varlog <- sdlog^2
            omega <- exp(varlog)
            H <- (omega * (1 + varlog) - (1 + 2 * varlog))^-1
            var.threshold <- (varlog/n) * (beta2/omega) * H
            if (ci.parameter == "threshold") {
                ci.obj <- ci.normal.approx(threshold, sqrt(var.threshold), 
                  n = n, df = n - 2, ci.type = ci.type, alpha = alpha)
                ci.obj$parameter <- "threshold"
            } else {
                var.beta <- (varlog/n) * beta2 * (1 + H)
                cov.threshold.beta <- -((sdlog^3)/n) * (beta2/sqrt(omega)) * 
                  H
                median <- threshold + beta
                var.median <- var.threshold + var.beta + 2 * 
                  cov.threshold.beta
                ci.obj <- ci.normal.approx(median, sqrt(var.median), 
                  n = n, df = n - 2, ci.type = ci.type, alpha = alpha)
                ci.obj$parameter <- "median"
            }
            ci.obj$method <- paste("Normal Approximation\n", 
                space(33), "Based on Asymptotic Variance", sep = "")
        }, likelihood.profile = {
            ci.obj <- ci.lnorm3.likelihood.profile(threshold = threshold, 
                meanlog = meanlog, sdlog = sdlog, x = x, x1 = x1, 
                n = n, eta.lb = eta.lb, eta.ub = eta.ub, ci.type = ci.type, 
                alpha = alpha, ci.parameter = ci.parameter)
        }, skewness = {
            ci.obj <- ci.lnorm3.zero.skew(threshold = threshold, 
                meanlog = meanlog, sdlog = sdlog, x = x, x1 = x1, 
                n = n, threshold.lb = threshold.lb, threshold.ub = threshold.ub, 
                ci.type = ci.type, alpha = alpha, ci.parameter = ci.parameter)
        })
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
