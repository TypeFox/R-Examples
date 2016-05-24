enormMultiplyCensored <-
function (x, censored, method = "mle", censoring.side = "left", 
    ci = FALSE, ci.method = "profile.likelihood", ci.type = "two-sided", 
    conf.level = 0.95, n.bootstraps = 1000, pivot.statistic = "z", 
    nmc = 1000, seed = NULL, ...) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!is.vector(censored, mode = "numeric") && !is.vector(censored, 
        mode = "logical")) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    data.name <- deparse(substitute(x))
    censoring.name <- deparse(substitute(censored))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'censored' removed."))
    }
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    n.cen <- sum(censored)
    if (n.cen == 0) 
        stop("No censored values indicated by 'censored'.")
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    N <- length(x)
    method <- match.arg(method, c("mle", "qq.reg", "impute.w.qq.reg", 
        "half.cen.level"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    x.cen <- x[censored]
    c.vec <- table(x.cen)
    cen.levels <- sort(unique(x.cen))
    K <- length(cen.levels)
    if (method == "half.cen.level" && (censoring.side == "right" || 
        any(cen.levels < 2 * .Machine$double.eps))) 
        stop(paste("The method 'half.cen.level' is applicable only for", 
            "left-censored data with positive censoring levels"))
    ci.method <- match.arg(ci.method, c("normal.approx", "bootstrap", 
        "profile.likelihood", "gpq"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    if (ci && ci.method == "profile.likelihood" && method != 
        "mle") 
        stop("When ci.method=\"profile.likelihood\" you must set method=\"mle\"")
    pivot.statistic <- match.arg(pivot.statistic, c("z", "t"))
    est.fcn <- paste("enormMultiplyCensored", method, sep = ".")
    if (!ci || !(ci.method %in% c("bootstrap", "gpq"))) {
        param.ci.list <- do.call(est.fcn, list(x = x, censored = censored, 
            N = N, cen.levels = cen.levels, K = K, c.vec = c.vec, 
            n.cen = n.cen, censoring.side = censoring.side, ci = ci, 
            ci.method = ci.method, ci.type = ci.type, conf.level = conf.level, 
            pivot.statistic = pivot.statistic, ...))
    }
    else {
        if (ci.method == "bootstrap") {
            param.list <- do.call(est.fcn, list(x = x, censored = censored, 
                N = N, cen.levels = cen.levels, K = K, c.vec = c.vec, 
                n.cen = n.cen, censoring.side = censoring.side, 
                ci = FALSE, ci.method = ci.method, ci.type = ci.type, 
                conf.level = conf.level, ...))
            ci.list <- enormMultiplyCensored.bootstrap.ci(x = x, 
                censored = censored, N = N, cen.levels = cen.levels, 
                K = K, c.vec = c.vec, n.cen = n.cen, censoring.side = censoring.side, 
                est.fcn = est.fcn, ci.type = ci.type, conf.level = conf.level, 
                n.bootstraps = n.bootstraps, obs.mean = param.list$parameters["mean"], 
                ...)
            param.ci.list <- c(param.list, list(ci.obj = ci.list))
        }
        else {
            diffs <- diff(sort(x))
            const <- min(diffs[diffs > 0])/2
            if (censoring.side == "right") 
                const <- -const
            new.x <- x
            new.x[!censored] <- new.x[!censored] + const
            new.censored <- censored[order(new.x)]
            cen.index <- (1:N)[new.censored]
            param.list <- do.call(est.fcn, list(x = x, censored = censored, 
                N = N, cen.levels = cen.levels, K = K, c.vec = c.vec, 
                n.cen = n.cen, censoring.side = censoring.side, 
                ci = FALSE, ci.method = ci.method, ci.type = ci.type, 
                conf.level = conf.level, ...))
            params <- param.list$parameters
            alpha <- 1 - conf.level
            probs <- switch(ci.type, `two-sided` = c((1 + conf.level)/2, 
                alpha/2), lower = conf.level, upper = alpha)
            gpq <- gpqCiNormMultiplyCensored(n = N, cen.index = cen.index, 
                probs = probs, nmc = nmc, method = method, censoring.side = censoring.side, 
                seed = seed, names = FALSE)
            limits <- switch(ci.type, `two-sided` = params["mean"] - 
                gpq * params["sd"], lower = c(params["mean"] - 
                gpq * params["sd"], Inf), upper = c(-Inf, params["mean"] - 
                gpq * params["sd"]))
            names(limits) <- c("LCL", "UCL")
            ci.list <- list(name = "Confidence", parameter = "mean", 
                limits = limits, type = ifelse(ci.type == "two.sided", 
                  "two-sided", ci.type), method = "Generalized Pivotal Quantity", 
                conf.level = conf.level, nmc = nmc)
            oldClass(ci.list) <- "intervalEstimateCensored"
            param.ci.list <- c(param.list, list(ci.obj = ci.list))
        }
    }
    method.string <- switch(method, mle = "MLE", qq.reg = "Q-Q Regression (ROS)", 
        impute.w.qq.reg = "Imputation With Q-Q Regression (ROS)", 
        half.cen.level = "Half Censoring Level")
    ret.list <- list(distribution = "Normal", sample.size = N, 
        censoring.side = censoring.side, censoring.levels = cen.levels, 
        percent.censored = (100 * n.cen)/N, parameters = param.ci.list$parameters, 
        n.param.est = 2, method = method.string)
    if (method %in% c("qq.reg", "impute.w.qq.reg")) 
        ret.list <- c(ret.list, list(prob.method = param.ci.list$prob.method, 
            plot.pos.con = param.ci.list$plot.pos.con))
    ret.list <- c(ret.list, list(data.name = data.name, censoring.name = censoring.name, 
        bad.obs = bad.obs))
    if (ci) {
        ret.list <- c(ret.list, list(interval = param.ci.list$ci.obj))
        if (!is.null(param.ci.list$var.cov.params)) 
            ret.list <- c(ret.list, list(var.cov.params = param.ci.list$var.cov.params))
    }
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
