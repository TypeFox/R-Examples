eqnormCensored <-
function (x, censored, censoring.side = "left", p = 0.5, method = "mle", 
    ci = FALSE, ci.method = "exact.for.complete", ci.type = "two-sided", 
    conf.level = 0.95, digits = 0, nmc = 1000, seed = NULL) 
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
    if (!is.vector(p, mode = "numeric")) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1")
    ci.method <- match.arg(ci.method, c("exact.for.complete", 
        "gpq", "normal.approx"))
    ret.list <- enormCensored(x, censored = censored, censoring.side = censoring.side, 
        method = method, ci = FALSE)
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    params <- ret.list$parameters
    q <- qnorm(p, mean = params["mean"], sd = params["sd"])
    if (length(p) == 1 && p == 0.5) 
        names(q) <- "Median"
    else {
        pct <- round(100 * p, digits)
        names(q) <- paste(pct, number.suffix(pct), " %ile", sep = "")
    }
    ret.list <- c(ret.list, list(quantiles = q))
    ret.list$quantile.method <- paste("Quantile(s) Based on\n", 
        space(33), ret.list$method, " Estimators", sep = "")
    if (ci) {
        if (length(p) > 1 || p <= 0 || p >= 1) 
            stop(paste("When 'ci' = TRUE, 'p' must be a scalar", 
                "larger than 0 and less than 1."))
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        if (!is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || !is.finite(conf.level) || conf.level <= 0 || 
            conf.level >= 1) 
            stop("'conf.level' must be a numeric scalar between 0 and 1")
        n <- length(x)
        if (ci.method != "gpq") {
            ci.method.arg <- ci.method
            if (ci.method == "exact.for.complete") 
                ci.method.arg <- "exact"
            ci.obj <- ci.qnorm(p = p, muhat = params["mean"], 
                sdhat = params["sd"], n = n, method = ci.method.arg, 
                ci.type = ci.type, alpha = 1 - conf.level, digits = digits)
            if (ci.method == "exact.for.complete" || (ci.method == 
                "normal.approx" && p == 0.5)) 
                ci.obj$method <- paste("Exact for\n", space(33), 
                  "Complete Data", sep = "")
        }
        else {
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
            probs <- switch(ci.type, lower = 1 - conf.level, 
                upper = conf.level, `two-sided` = c((1 - conf.level)/2, 
                  (1 + conf.level)/2))
            if (multiple) {
                diffs <- diff(sort(x))
                const <- min(diffs[diffs > 0])/2
                if (censoring.side == "right") 
                  const <- -const
                new.x <- x
                new.x[!censored] <- new.x[!censored] + const
                new.censored <- censored[order(new.x)]
                cen.index <- (1:n)[new.censored]
                if (ci.type == "lower") {
                  gpq <- gpqTolIntNormMultiplyCensored(n = n, 
                    cen.index = cen.index, p = p, probs = probs, 
                    nmc = nmc, censoring.side = censoring.side, 
                    seed = seed, names = FALSE)
                  limits <- c(params["mean"] + gpq * params["sd"], 
                    Inf)
                }
                else if (ci.type == "upper") {
                  gpq <- gpqTolIntNormMultiplyCensored(n = n, 
                    cen.index = cen.index, p = p, probs = probs, 
                    nmc = nmc, censoring.side = censoring.side, 
                    seed = seed, names = FALSE)
                  limits <- c(-Inf, params["mean"] + gpq * params["sd"])
                }
                else {
                  gpq.lower <- gpqTolIntNormMultiplyCensored(n = n, 
                    cen.index = cen.index, p = p, probs = probs[1], 
                    nmc = nmc, censoring.side = censoring.side, 
                    seed = seed, names = FALSE)
                  gpq.upper <- gpqTolIntNormMultiplyCensored(n = n, 
                    cen.index = cen.index, p = p, probs = probs[2], 
                    nmc = nmc, censoring.side = censoring.side, 
                    seed = seed, names = FALSE)
                  limits <- c(params["mean"] + gpq.lower * params["sd"], 
                    params["mean"] + gpq.upper * params["sd"])
                }
            }
            else {
                if (ci.type == "lower") {
                  gpq <- gpqTolIntNormSinglyCensored(n = n, n.cen = n.cen, 
                    p = p, probs = probs, nmc = nmc, censoring.side = censoring.side, 
                    seed = seed, names = FALSE)
                  limits <- c(params["mean"] + gpq * params["sd"], 
                    Inf)
                }
                else if (ci.type == "upper") {
                  gpq <- gpqTolIntNormSinglyCensored(n = n, n.cen = n.cen, 
                    p = p, probs = probs, nmc = nmc, censoring.side = censoring.side, 
                    seed = seed, names = FALSE)
                  limits <- c(-Inf, params["mean"] + gpq * params["sd"])
                }
                else {
                  gpq.lower <- gpqTolIntNormSinglyCensored(n = n, 
                    n.cen = n.cen, p = p, probs = probs[1], nmc = nmc, 
                    censoring.side = censoring.side, seed = seed, 
                    names = FALSE)
                  gpq.upper <- gpqTolIntNormSinglyCensored(n = n, 
                    n.cen = n.cen, p = p, probs = probs[2], nmc = nmc, 
                    censoring.side = censoring.side, seed = seed, 
                    names = FALSE)
                  limits <- c(params["mean"] + gpq.lower * params["sd"], 
                    params["mean"] + gpq.upper * params["sd"])
                }
            }
            names(limits) <- c("LCL", "UCL")
            ci.obj <- list(name = "Confidence", parameter = names(q), 
                limits = limits, type = ifelse(ci.type == "two.sided", 
                  "two-sided", ci.type), method = "Generalized Pivotal Quantity", 
                conf.level = conf.level, nmc = nmc)
        }
        oldClass(ci.obj) <- "intervalEstimateCensored"
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
