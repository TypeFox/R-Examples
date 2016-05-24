elnormAlt <-
function (x, method = "mvue", ci = FALSE, ci.type = "two-sided", 
    ci.method = "land", conf.level = 0.95, parkin.list = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2 || any(x <= 0)) 
        stop(paste("'x' must contain at least 2 non-missing distinct values,", 
            "and all non-missing values must be positive."))
    method <- match.arg(method, c("mvue", "qmle", "mle", "mme", 
        "mmue"))
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method, c("land", "zou", "parkin", 
            "cox", "normal.approx"))
        if (ci.method == "land" && n < 3) 
            stop(paste("When ci=TRUE and ci.method=\"land\",", 
                "'x' must contain at least 3 values"))
    }
    meanlog <- mean(log(x))
    s2 <- var(log(x))
    sdlog <- sqrt(s2)
    s2.mle <- ((n - 1)/n) * s2
    df <- n - 1
    switch(method, mvue = {
        muhat <- exp(meanlog) * finneys.g(n - 1, s2/2)
        sdhat <- sqrt(exp(2 * meanlog) * (finneys.g(n - 1, 2 * 
            s2) - finneys.g(n - 1, (s2 * (n - 2))/(n - 1))))
        if (ci && ci.method == "normal.approx") sd.muhat <- sqrt(exp(2 * 
            meanlog) * ((finneys.g(n - 1, s2/2)^2) - finneys.g(n - 
            1, (s2 * (n - 2))/(n - 1))))
    }, qmle = {
        muhat <- exp(meanlog + s2/2)
        sdhat <- muhat * sqrt(exp(s2) - 1)
        if (ci && ci.method == "normal.approx") sd.muhat <- sqrt(exp(2 * 
            meanlog + s2/n) * (exp(s2/n) * ((1 - (2 * s2)/df)^(-df/2)) - 
            ((1 - s2/df)^(-df))))
    }, mle = {
        muhat <- exp(meanlog + s2.mle/2)
        sdhat <- muhat * sqrt(exp(s2.mle) - 1)
        if (ci && ci.method == "normal.approx") sd.muhat <- sqrt(exp(2 * 
            meanlog + s2/n) * (exp(s2/n) * ((1 - (2 * s2)/n)^(-df/2)) - 
            ((1 - s2/n)^(-df))))
    }, mme = , mmue = {
        muhat <- mean(x)
        sdhat <- ifelse(method != "mme", sd(x), sqrt((n - 1)/n) * 
            sd(x))
        if (ci && ci.method == "normal.approx") sd.muhat <- sdhat/sqrt(n)
    })
    dist.params <- c(mean = muhat, cv = sdhat/muhat)
    ret.list <- list(distribution = "Lognormal", sample.size = n, 
        parameters = dist.params, n.param.est = 2, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        switch(ci.method, land = {
            ci.obj <- ci.lnorm.land(meanlog, sdlog, n, ci.type, 
                alpha = 1 - conf.level)
        }, zou = {
            ci.obj <- ci.lnorm.zou(meanlog, sdlog, n, ci.type, 
                alpha = 1 - conf.level)
        }, parkin = {
            p.hat <- pnorm(sdlog/2)
            if (is.null(parkin.list)) {
                parkin.list <- list(lcl.rank = NULL, ucl.rank = NULL, 
                  ci.method = ifelse(n <= 20, "exact", "normal.approx"), 
                  approx.conf.level = conf.level)
            }
            ci.obj <- do.call("eqnpar", c(list(x = x, p = p.hat, 
                ci = TRUE, lb = 0, ci.type = ci.type), parkin.list))$interval
            ci.obj$parameter <- "mean"
            ci.obj$method <- "Parkin"
        }, cox = {
            beta.hat <- meanlog + (s2/2)
            sd.beta.hat <- sqrt(s2/n + (s2^2)/(2 * (n + 1)))
            ci.obj <- ci.normal.approx(beta.hat, sd.beta.hat, 
                n, df, ci.type, alpha = 1 - conf.level)
            ci.obj$limits <- exp(ci.obj$limits)
            ci.obj$parameter <- "mean"
            ci.obj$method <- "Cox"
        }, normal.approx = {
            ci.obj <- ci.normal.approx(muhat, sd.muhat, n, df, 
                ci.type, alpha = 1 - conf.level, lb = 0)
            ci.obj$parameter <- "mean"
        })
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
