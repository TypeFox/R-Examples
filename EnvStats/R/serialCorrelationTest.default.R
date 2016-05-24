serialCorrelationTest.default <-
function (x, test = "rank.von.Neumann", alternative = "two.sided", 
    conf.level = 0.95, ...) 
{
    data.name <- deparse(substitute(x))
    test <- match.arg(test, c("rank.von.Neumann", "AR1.yw", "AR1.mle"))
    is.ts <- is.ts(x)
    if (is.ts) {
        ts.x <- x
        x <- unclass(x)
        if (!is.null(dim(x))) 
            stop("When 'x' is an object of class 'ts', it must be a univariate time series")
        attributes(x) <- NULL
        if (!is.vector(x, mode = "numeric")) 
            stop("When 'x' is an object of class 'ts', it must be a univariate numeric time series")
    }
    else {
        if (!is.vector(x, mode = "numeric")) 
            stop("'x' must either be a vector of class 'ts' or a numeric vector.")
    }
    n <- length(x)
    if (test == "AR1.mle") {
        if (any(is.infinite(x)) || any(is.nan(x))) 
            stop("'x' cannot contain infinite (Inf, -Inf) or undefined (NaN) values")
        bad.obs <- sum(is.na(x))
        if ((n - bad.obs) < 3) 
            stop("'x' must have at least 3 non-missing values.")
    }
    else {
        if (!all(is.finite(x))) 
            stop(paste("'x' cannot contain missing (NA),", "infinite (Inf, -Inf), or undefined (NaN) values."))
        bad.obs <- 0
        if (n < 3) 
            stop("'x' must have at least 3 values.")
    }
    alternative <- match.arg(alternative, c("two.sided", "greater", 
        "less"))
    if (test != "AR1.mle") {
        estimate <- ar.yw(x, aic = FALSE, order = 1)$ar
        names(estimate) <- "rho"
        var.estimate <- (1 - estimate^2)/n
        estimation.method <- "Yule-Walker"
    }
    sep.string <- paste("\n", space(33), sep = "")
    switch(test, rank.von.Neumann = {
        method <- paste("Rank von Neumann Test for", "Lag-1 Autocorrelation", 
            sep = sep.string)
        ties <- rle(sort(x))$lengths
        pct.ties <- (100 * sum(ties[ties > 1]))/n
        if (pct.ties > 0) warning(paste(round(pct.ties, 0), "% of the observations are tied.  ", 
            "The p-value may not be accurate.", sep = ""))
        rvn <- 12/(n * (n^2 - 1)) * sum(diff(rank(x))^2)
        if (n > 100) {
            z <- (rvn - 2)/sqrt(20/(5 * n + 7))
            p.value <- switch(alternative, two.sided = 2 * (1 - 
                pnorm(abs(z))), greater = pnorm(z), less = 1 - 
                pnorm(z))
            statistic <- z
            names(statistic) <- "z"
            method <- paste(method, "(Normal Approximation)", 
                sep = sep.string)
        } else if (n > 10) {
            shape1 <- ((5 * n * (n + 1) * (n - 1)^2)/(2 * (n - 
                2) * (5 * n^2 - 2 * n - 9))) - 0.5
            p.value <- switch(alternative, two.sided = {
                if (rvn <= 2) {
                  pbeta(rvn/4, shape1, shape1) + 1 - pbeta((4 - 
                    rvn)/4, shape1, shape1)
                } else {
                  pbeta((4 - rvn)/4, shape1, shape1) + 1 - pbeta(rvn/4, 
                    shape1, shape1)
                }
            }, greater = pbeta(rvn/4, shape1, shape1), less = 1 - 
                pbeta(rvn/4, shape1, shape1))
            statistic <- rvn
            names(statistic) <- "RVN"
            method <- paste(method, "(Beta Approximation)", sep = sep.string)
        } else {
            dist.mat <- get(paste(".rank.von.neumann.dist.mat.", 
                n, sep = ""))
            ord.stats <- dist.mat[, "Order.Statistics"]
            dens <- dist.mat[, "Densities"]
            cum.probs <- dist.mat[, "Cumulative.Probabilities"]
            new.rvn <- approx(ord.stats, ord.stats, rvn, method = "constant", 
                rule = 2, f = 1)$y
            n.ord.stats <- length(ord.stats)
            index <- (1:n.ord.stats)[abs(new.rvn - ord.stats) < 
                1e-06]
            p.value <- switch(alternative, two.sided = {
                if (index <= (n.ord.stats + 1)/2) {
                  cum.probs[index] + 1 - cum.probs[n.ord.stats - 
                    index + 1] + dens[n.ord.stats - index + 1]
                } else {
                  cum.probs[n.ord.stats - index + 1] + 1 - cum.probs[index] + 
                    dens[index]
                }
            }, greater = cum.probs[index], less = 1 - cum.probs[index] + 
                dens[index])
            statistic <- rvn
            names(statistic) <- "RVN"
            method <- paste(method, "(Exact Method)", sep = sep.string)
        }
    }, AR1.yw = {
        z <- sqrt(n) * estimate
        p.value <- switch(alternative, two.sided = 2 * (1 - pnorm(abs(z))), 
            greater = 1 - pnorm(z), less = pnorm(z))
        statistic <- z
        names(statistic) <- "z"
        method <- paste("z-Test for", "Lag-1 Autocorrelation", 
            "(Normal Approximation)", sep = sep.string)
    }, AR1.mle = {
        if (is.ts) x <- ts.x
        ar1.list <- arima(x, order = c(1, 0, 0), include.mean = TRUE, 
            method = "ML")
        estimate <- ar1.list$coef["ar1"]
        var.estimate <- ar1.list$var.coef["ar1", "ar1"]
        names(estimate) <- "rho"
        estimation.method <- "Maximum Likelihood"
        z <- estimate/sqrt(var.estimate)
        p.value <- switch(alternative, two.sided = 2 * (1 - pnorm(abs(z))), 
            greater = 1 - pnorm(z), less = pnorm(z))
        statistic <- z
        names(statistic) <- "z"
        method <- paste("z-Test for", "Lag-1 Autocorrelation", 
            "(Wald Test Based on MLE)", sep = sep.string)
    })
    null.value <- 0
    names(null.value) <- "rho"
    ret.list <- list(statistic = statistic, parameters = NULL, 
        p.value = p.value, estimate = estimate, estimation.method = estimation.method, 
        null.value = null.value, alternative = alternative, method = method, 
        sample.size = n, data.name = data.name, bad.obs = bad.obs)
    ci.type <- switch(alternative, two.sided = "two.sided", greater = "lower", 
        less = "upper")
    ci.obj <- ci.normal.approx(theta.hat = estimate, sd.theta.hat = sqrt(var.estimate), 
        n = n, ci.type = ci.type, alpha = 1 - conf.level, lb = -1, 
        ub = 1, test.statistic = "z")
    ci.obj$parameter <- "rho"
    ret.list <- c(ret.list, list(interval = ci.obj))
    oldClass(ret.list) <- "htest"
    ret.list
}
