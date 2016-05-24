enormSinglyCensored <-
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
    method <- match.arg(method, c("mle", "bcmle", "qq.reg", "qq.reg.w.cen.level", 
        "impute.w.qq.reg", "impute.w.qq.reg.w.cen.level", "impute.w.mle", 
        "iterative.impute.w.qq.reg", "m.est", "half.cen.level"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    T1 <- unique(x[censored])
    if (length(T1) > 1) 
        stop(paste("More than one censoring level.  Use 'enormMultiplyCensored'."))
    if (censoring.side == "left") {
        if (T1 > min(x.no.cen)) 
            stop(paste("For singly left-censored data,", "all uncensored observations must be bigger than", 
                "or equal to the censoring level. ", "Use enormMultiplyCensored."))
    }
    else {
        if (T1 < max(x.no.cen)) 
            stop(paste("For singly right-censored data,", "all uncensored observations must be less than", 
                "or equal to the censoring level. ", "Use enormMultiplyCensored."))
    }
    if (method == "half.cen.level" && (censoring.side == "right" || 
        T1 < 2 * .Machine$double.eps)) 
        stop(paste("The method 'half.cen.level' is applicable only for", 
            "left-censored data with a positive censoring level"))
    ci.method <- match.arg(ci.method, c("normal.approx", "normal.approx.w.cov", 
        "bootstrap", "profile.likelihood", "gpq"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    if (ci && ci.method == "profile.likelihood" && method != 
        "mle") 
        stop("When ci.method=\"profile.likelihood\" you must set method=\"mle\"")
    if (ci && ci.method == "normal.approx.w.cov" && !(method %in% 
        c("mle", "bcmle"))) 
        stop(paste("When ci.method='normal.approx.w.cov'", "you must set method='mle' or method='bcmle'"))
    est.fcn <- paste("enormSinglyCensored", method, sep = ".")
    if (!ci || !(ci.method %in% c("bootstrap", "gpq"))) {
        pivot.statistic <- match.arg(pivot.statistic, c("z", 
            "t"))
        param.ci.list <- do.call(est.fcn, list(x = x, censored = censored, 
            N = N, T1 = T1, n.cen = n.cen, censoring.side = censoring.side, 
            ci = ci, ci.method = ci.method, ci.type = ci.type, 
            conf.level = conf.level, pivot.statistic = pivot.statistic, 
            ...))
    }
    else {
        if (ci.method == "bootstrap") {
            param.list <- do.call(est.fcn, list(x = x, censored = censored, 
                N = N, T1 = T1, n.cen = n.cen, censoring.side = censoring.side, 
                ci = FALSE, ci.method = ci.method, ci.type = ci.type, 
                conf.level = conf.level, ...))
            ci.list <- enormSinglyCensored.bootstrap.ci(x = x, 
                censored = censored, N = N, T1 = T1, censoring.side = censoring.side, 
                est.fcn = est.fcn, ci.type = ci.type, conf.level = conf.level, 
                n.bootstraps = n.bootstraps, obs.mean = param.list$parameters["mean"], 
                ...)
            param.ci.list <- c(param.list, list(ci.obj = ci.list))
        }
        else {
            param.list <- do.call(est.fcn, list(x = x, censored = censored, 
                N = N, T1 = T1, n.cen = n.cen, censoring.side = censoring.side, 
                ci = FALSE, ci.method = ci.method, ci.type = ci.type, 
                conf.level = conf.level, ...))
            params <- param.list$parameters
            alpha <- 1 - conf.level
            probs <- switch(ci.type, `two-sided` = c((1 + conf.level)/2, 
                alpha/2), lower = conf.level, upper = alpha)
            gpq <- gpqCiNormSinglyCensored(n = N, n.cen = n.cen, 
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
    if (method == "m.est") {
        arg.list <- list(...)
        if (is.null(arg.list) || !is.element("t.df", names(arg.list))) 
            t.df <- 3
        method <- paste("m.est (t.df=", t.df, ")", sep = "")
    }
    method.string <- switch(method, mle = "MLE", bcmle = "Bias-corrected MLE", 
        qq.reg = "Q-Q Regression (ROS)", qq.reg.w.cen.level = paste("Q-Q Regression (ROS)\n", 
            space(33), "with Censoring Level", sep = ""), impute.w.qq.reg = paste("Imputation with\n", 
            space(33), "Q-Q Regression (ROS)", sep = ""), impute.w.qq.reg.w.cen.level = paste("Imputation with\n", 
            space(33), "Q-Q Regression (ROS)\n", space(33), "with Censoring Level", 
            sep = ""), impute.w.mle = "Imputation with MLE", 
        iterative.impute.w.qq.reg = paste("Iterative Imputation\n", 
            space(33), "with Q-Q Regression (ROS)", sep = ""), 
        half.cen.level = "Half Censoring Level", m.est = "M-Estimators")
    ret.list <- list(distribution = "Normal", sample.size = N, 
        censoring.side = censoring.side, censoring.levels = T1, 
        percent.censored = (100 * n.cen)/N, parameters = param.ci.list$parameters, 
        n.param.est = 2, method = method.string, data.name = data.name, 
        censoring.name = censoring.name, bad.obs = bad.obs)
    if (ci) {
        ret.list <- c(ret.list, list(interval = param.ci.list$ci.obj))
        if (!is.null(param.ci.list$var.cov.params)) 
            ret.list <- c(ret.list, list(var.cov.params = param.ci.list$var.cov.params))
    }
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
