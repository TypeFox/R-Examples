enparCensored <-
function (x, censored, censoring.side = "left", correct.se = FALSE, 
    left.censored.min = "DL", right.censored.max = "DL", ci = FALSE, 
    ci.method = "normal.approx", ci.type = "two-sided", conf.level = 0.95, 
    pivot.statistic = "z", ci.sample.size = NULL, n.bootstraps = 1000) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!is.vector(censored, mode = "numeric") & !is.vector(censored, 
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
    if (any(x <= 0)) 
        stop("All values of 'x' must be positive")
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
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    x.cen <- x[censored]
    cen.levels <- sort(unique(x.cen))
    ci.method <- match.arg(ci.method, c("normal.approx", "bootstrap"))
    ci.type <- match.arg(ci.type, c("two-sided", "lower", "upper"))
    pivot.statistic <- match.arg(pivot.statistic, c("z", "t"))
    if (censoring.side == "left") {
        dl <- min(cen.levels)
        if (dl <= min(x.no.cen)) {
            if (length(left.censored.min) != 1) 
                stop("The argument 'left.censored.min' must have length 1")
            if (is.character(left.censored.min)) {
                if (!(left.censored.min %in% c("DL/2", "DL", 
                  "Ignore"))) 
                  stop(paste("When the argument 'left.censored.min' is a", 
                    "character string, it must equal 'DL', 'DL/2', or 'Ignore'"))
            }
            else {
                if (!is.numeric(left.censored.min) || left.censored.min > 
                  dl || left.censored.min <= 0) 
                  stop(paste("When 'left.censored.min' is a numeric scalar", 
                    "it must be less than or equal to the smallest", 
                    "censoring limit and greater than 0"))
            }
        }
    }
    else {
        dl <- max(cen.levels)
        if (dl >= max(x.no.cen)) {
            if (length(right.censored.max) != 1) 
                stop("The argument 'right.censored.max' must have length 1")
            if (is.character(right.censored.max)) {
                if (!(right.censored.max %in% c("Ignore", "DL"))) 
                  stop(paste("When the argument 'right.censored.max'", 
                    "is a character string, it must equal 'DL' or 'Ignore'"))
            }
            else {
                if (!is.numeric(right.censored.max) || right.censored.max < 
                  dl) 
                  stop(paste("When 'right.censored.max' is a numeric scalar", 
                    "it must be greater than or equal to the largest", 
                    "censoring limit"))
            }
        }
    }
    if (!ci || ci.method != "bootstrap") {
        param.ci.list <- enparCensored.km(x = x, censored = censored, 
            censoring.side = censoring.side, correct.se = correct.se, 
            left.censored.min = left.censored.min, right.censored.max = right.censored.max, 
            ci = ci, ci.type = ci.type, conf.level = conf.level, 
            pivot.statistic = pivot.statistic, ci.sample.size = ci.sample.size)
    }
    else {
        param.ci.list <- enparCensored.km(x = x, censored = censored, 
            censoring.side = censoring.side, correct.se = correct.se, 
            left.censored.min = left.censored.min, right.censored.max = right.censored.max, 
            ci = FALSE)
        ci.list <- enparCensored.bootstrap.ci(x = x, censored = censored, 
            censoring.side = censoring.side, correct.se = correct.se, 
            left.censored.min = left.censored.min, right.censored.max = right.censored.max, 
            est.fcn = "enparCensored.km", ci.type = ci.type, 
            conf.level = conf.level, n.bootstraps = n.bootstraps, 
            obs.mean = param.ci.list$parameters["mean"], obs.se.mean = param.ci.list$parameters["se.mean"])
        param.ci.list <- c(param.ci.list, list(ci.obj = ci.list))
    }
    method <- "Kaplan-Meier"
    if (correct.se) 
        method <- paste(method, "\n", space(33), "(Bias-corrected se.mean)", 
            sep = "")
    ret.list <- list(distribution = "None", sample.size = N, 
        censoring.side = censoring.side, censoring.levels = cen.levels, 
        percent.censored = (100 * n.cen)/N, parameters = param.ci.list$parameters, 
        n.param.est = 2, method = method, data.name = data.name, 
        censoring.name = censoring.name, bad.obs = bad.obs)
    if (ci) {
        ret.list <- c(ret.list, list(interval = param.ci.list$ci.obj))
        if (!is.null(param.ci.list$var.cov.params)) 
            ret.list <- c(ret.list, list(var.cov.params = param.ci.list$var.cov.params))
    }
    oldClass(ret.list) <- "estimateCensored"
    ret.list
}
