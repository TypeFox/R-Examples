tolIntNormCensored <-
function (x, censored, censoring.side = "left", coverage = 0.95, 
    cov.type = "content", ti.type = "two-sided", conf.level = 0.95, 
    method = "mle", ti.method = "exact.for.complete", seed = NULL, 
    nmc = 1000) 
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
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (!is.numeric(coverage) || length(coverage) != 1 || is.na(coverage) || 
        coverage <= 0 || coverage >= 1) 
        stop("'coverage' must be a scalar greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    ti.type <- match.arg(ti.type, c("two-sided", "lower", "upper"))
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    multiple <- TRUE
    T.vec <- unique(x[censored])
    if (length(T.vec) == 1) {
        if (censoring.side == "left") {
            if (T.vec <= min(x.no.cen)) 
                multiple <- FALSE
        }
        else {
            if (T.vec >= max(x.no.cen)) 
                multiple <- FALSE
        }
    }
    if (multiple) {
        method <- match.arg(method, c("mle", "qq.reg", "impute.w.qq.reg", 
            "half.cen.level"))
    }
    else {
        method <- match.arg(method, c("mle", "bcmle", "qq.reg", 
            "qq.reg.w.cen.level", "impute.w.qq.reg", "impute.w.qq.reg.w.cen.level", 
            "impute.w.mle", "iterative.impute.w.qq.reg", "m.est", 
            "half.cen.level"))
    }
    ti.method <- match.arg(ti.method, c("exact.for.complete", 
        "wald.wolfowitz.for.complete", "gpq"))
    if (ti.method == "gpq" && cov.type != "content") 
        stop("When ti.method='gpq' you must set cov.type='content'")
    ret.list <- enormCensored(x, censored = censored, method = method, 
        censoring.side = censoring.side, ci = FALSE)
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    if (ti.method != "gpq") {
        ti.method.arg <- switch(ti.method, exact.for.complete = "exact", 
            wald.wolfowitz.for.complete = "wald.wolfowitz")
        ret.list <- tolIntNorm(x = ret.list, coverage = coverage, 
            cov.type = cov.type, ti.type = ti.type, conf.level = conf.level, 
            method = ti.method.arg)
        ret.list$interval$method <- paste(ret.list$interval$method, 
            " for\n", space(33), "Complete Data", sep = "")
        oldClass(ret.list$interval) <- "intervalEstimateCensored"
    }
    else {
        n <- length(x)
        params <- ret.list$parameters
        probs <- switch(ti.type, lower = 1 - conf.level, upper = conf.level, 
            `two-sided` = c((1 - conf.level)/2, (1 + conf.level)/2))
        p <- switch(ti.type, lower = 1 - coverage, upper = coverage, 
            `two-sided` = c((1 - coverage)/2, (1 + coverage)/2))
        if (multiple) {
            diffs <- diff(sort(x))
            const <- min(diffs[diffs > 0])/2
            if (censoring.side == "right") 
                const <- -const
            new.x <- x
            new.x[!censored] <- new.x[!censored] + const
            new.censored <- censored[order(new.x)]
            cen.index <- (1:n)[new.censored]
            if (ti.type == "lower") {
                gpq <- gpqTolIntNormMultiplyCensored(n = n, cen.index = cen.index, 
                  p = p, probs = probs, nmc = nmc, method = method, 
                  censoring.side = censoring.side, seed = seed, 
                  names = FALSE)
                limits <- c(params["mean"] + gpq * params["sd"], 
                  Inf)
            }
            else if (ti.type == "upper") {
                gpq <- gpqTolIntNormMultiplyCensored(n = n, cen.index = cen.index, 
                  p = p, probs = probs, nmc = nmc, method = method, 
                  censoring.side = censoring.side, seed = seed, 
                  names = FALSE)
                limits <- c(-Inf, params["mean"] + gpq * params["sd"])
            }
            else {
                gpq.lower <- gpqTolIntNormMultiplyCensored(n = n, 
                  cen.index = cen.index, p = p[1], probs = probs[1], 
                  nmc = nmc, method = method, censoring.side = censoring.side, 
                  seed = seed, names = FALSE)
                gpq.upper <- gpqTolIntNormMultiplyCensored(n = n, 
                  cen.index = cen.index, p = p[2], probs = probs[2], 
                  nmc = nmc, method = method, censoring.side = censoring.side, 
                  seed = seed, names = FALSE)
                limits <- c(params["mean"] + gpq.lower * params["sd"], 
                  params["mean"] + gpq.upper * params["sd"])
            }
        }
        else {
            if (ti.type == "lower") {
                gpq <- gpqTolIntNormSinglyCensored(n = n, n.cen = n.cen, 
                  p = p, probs = probs, nmc = nmc, method = method, 
                  censoring.side = censoring.side, seed = seed, 
                  names = FALSE)
                limits <- c(params["mean"] + gpq * params["sd"], 
                  Inf)
            }
            else if (ti.type == "upper") {
                gpq <- gpqTolIntNormSinglyCensored(n = n, n.cen = n.cen, 
                  p = p, probs = probs, nmc = nmc, method = method, 
                  censoring.side = censoring.side, seed = seed, 
                  names = FALSE)
                limits <- c(-Inf, params["mean"] + gpq * params["sd"])
            }
            else {
                gpq.lower <- gpqTolIntNormSinglyCensored(n = n, 
                  n.cen = n.cen, p = p[1], probs = probs[1], 
                  nmc = nmc, method = method, censoring.side = censoring.side, 
                  seed = seed, names = FALSE)
                gpq.upper <- gpqTolIntNormSinglyCensored(n = n, 
                  n.cen = n.cen, p = p[2], probs = probs[2], 
                  nmc = nmc, method = method, censoring.side = censoring.side, 
                  seed = seed, names = FALSE)
                limits <- c(params["mean"] + gpq.lower * params["sd"], 
                  params["mean"] + gpq.upper * params["sd"])
            }
        }
        names(limits) <- c("LTL", "UTL")
        ti.obj <- list(name = "Tolerance", coverage = coverage, 
            coverage.type = cov.type, limits = limits, type = ti.type, 
            method = "Generalized Pivotal Quantity", conf.level = conf.level, 
            nmc = nmc)
        oldClass(ti.obj) <- "intervalEstimateCensored"
        ret.list <- c(ret.list, list(interval = ti.obj))
    }
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
