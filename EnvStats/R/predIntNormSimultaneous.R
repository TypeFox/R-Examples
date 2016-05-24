predIntNormSimultaneous <-
function (x, n.mean = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    delta.over.sigma = 0, pi.type = "upper", conf.level = 0.95, 
    K.tol = .Machine$double.eps^0.5) 
{
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"))
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    switch(rule, k.of.m = {
        if (!is.vector(k, mode = "numeric") || length(k) != 1 || 
            k != trunc(k) || k < 1 || !is.vector(n.mean, mode = "numeric") || 
            length(n.mean) != 1 || n.mean != trunc(n.mean) || 
            n.mean < 1 || !is.vector(m, mode = "numeric") || 
            length(m) != 1 || m != trunc(m) || m < 1 || !is.vector(r, 
            mode = "numeric") || length(r) != 1 || r != trunc(r) || 
            r < 1 || k > m) stop(paste("'k', 'n.mean', 'm', and 'r' must be positive integers,", 
            "and 'k' must be between 1 and 'm'"))
    }, CA = {
        if (!is.vector(n.mean, mode = "numeric") || length(n.mean) != 
            1 || n.mean != trunc(n.mean) || n.mean < 1 || !is.vector(m, 
            mode = "numeric") || length(m) != 1 || m != trunc(m) || 
            m < 1 || !is.vector(r, mode = "numeric") || length(r) != 
            1 || r != trunc(r) || r < 1) stop("'n.mean', 'm', and 'r' must be positive integers")
    }, Modified.CA = {
        if (!is.vector(n.mean, mode = "numeric") || length(n.mean) != 
            1 || n.mean != trunc(n.mean) || n.mean < 1 || !is.vector(m, 
            mode = "numeric") || length(m) != 1 || m != trunc(m) || 
            m < 1 || !is.vector(r, mode = "numeric") || length(r) != 
            1 || r != trunc(r) || r < 1) stop("'n.mean', 'm', and 'r' must be positive integers")
        m <- 4
    })
    if (!is.vector(delta.over.sigma, mode = "numeric") || length(delta.over.sigma) != 
        1 || !is.finite(delta.over.sigma)) 
        stop("'delta.over.sigma' must be a finite numeric scalar.")
    if (!is.vector(conf.level, mode = "numeric") || length(conf.level) != 
        1 || conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored") {
        if (x$distribution != "Normal") 
            stop(paste("'predIntNormSimultaneous' creates prediction intervals", 
                "for a normal distribution.  You have supplied an object", 
                "that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        xbar <- x$parameters["mean"]
        s <- x$parameters["sd"]
        n <- x$sample.size
        ret.list <- x
    }
    else {
        if (!is.vector(x, mode = "numeric")) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric vector"))
        data.name <- deparse(substitute(x))
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        n <- length(x)
        if (n < 2 || length(unique(x)) < 2) 
            stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
                "This is not true for 'x' =", data.name))
        ret.list <- enorm(x)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        xbar <- ret.list$parameters["mean"]
        s <- ret.list$parameters["sd"]
    }
    df <- n - 1
    K <- predIntNormSimultaneousK(n = n, df = df, n.mean = n.mean, 
        k = k, m = m, r = r, rule = rule, delta.over.sigma = delta.over.sigma, 
        pi.type = pi.type, conf.level = conf.level, K.tol = K.tol)
    switch(pi.type, two.sided = {
        LPL <- xbar - K * s
        UPL <- xbar + K * s
    }, lower = {
        LPL <- xbar - K * s
        UPL <- Inf
    }, upper = {
        LPL <- -Inf
        UPL <- xbar + K * s
    })
    limits <- c(LPL, UPL)
    names(limits) <- c("LPL", "UPL")
    string <- ifelse(rule == "k.of.m", "", paste(" (", rule, 
        " Rule)", sep = ""))
    if (rule == "CA") 
        k <- "First.or.all.of.next.m.minus.one"
    else if (rule == "Modified.CA") 
        k <- "First.or.at.least.two.of.next.three"
    pi.obj <- list(name = "Prediction", rule = rule, limits = limits, 
        type = ifelse(pi.type == "two.sided", "two-sided", pi.type), 
        method = paste("exact", string, sep = ""), conf.level = conf.level, 
        sample.size = n, dof = df, k = k, m = m, r = r, delta.over.sigma = delta.over.sigma, 
        n.mean = n.mean)
    oldClass(pi.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = pi.obj))
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    ret.list
}
