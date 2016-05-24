boxcoxSinglyCensored <-
function (x, censored, censoring.side = c("left", "right"), lambda = {
    if (optimize) 
        c(-2, 2)
    else seq(-2, 2, by = 0.5)
}, optimize = FALSE, objective.name = c("PPCC", "Shapiro-Wilk", 
    "Log-Likelihood"), eps = .Machine$double.eps, include.x.and.censored = TRUE) 
{
    data.name <- deparse(substitute(x))
    censoring.name <- deparse(substitute(censored))
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!((is.vector(censored, mode = "numeric") && !is.factor(censored)) || 
        is.vector(censored, mode = "logical"))) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    if (any(is.na(censored))) 
        stop("'censored' cannot contain missing values")
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    censoring.side <- match.arg(censoring.side)
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        censored <- censored[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    if (any(x <= 0)) 
        stop("All non-missing, finite values of 'x' must be positive")
    objective.name <- match.arg(objective.name)
    n.cen <- sum(censored)
    if (n.cen == 0) {
        warning(paste("No censored values indicated by 'censored',", 
            "so the function 'boxcox' was called."))
        ret.list <- boxcox(x = x, lambda = lambda, optimize = optimize, 
            objective.name = objective.name, eps = eps, include.x = include.x.and.censored)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        return(ret.list)
    }
    T1 <- unique(x[censored])
    if (length(T1) > 1) 
        stop(paste("More than one censoring level. ", "Use 'boxcoxMultiplyCensored'."))
    x.no.cen <- x[!censored]
    if (censoring.side == "left") {
        if (T1 > min(x.no.cen)) 
            stop(paste("For singly left-censored data,", "all uncensored observations must be bigger than", 
                "or equal to the censoring level. ", "Use 'boxcoxMultiplyCensored'."))
    }
    else {
        if (T1 < max(x.no.cen)) 
            stop(paste("For singly right-censored data,", "all uncensored observations must be less than", 
                "or equal to the censoring level. ", "Use 'boxcoxMultiplyCensored'."))
    }
    N <- length(x)
    if (!is.vector(lambda, mode = "numeric") || any(!is.finite(lambda))) 
        stop("'lambda' must be a numeric vector with no missing or infinite values")
    if (optimize && length(unique(lambda)) != 2) 
        stop(paste("When optimize=TRUE, 'lambda' must be a vector", 
            "with two unique values that specify the lower and", 
            "upper bounds for the optimization"))
    if (optimize && (1 < min(lambda) || 1 > max(lambda))) 
        stop("When optimize=TRUE, the range of 'lambda' must contain 1")
    lambda <- sort(lambda)
    objective.fcn <- switch(objective.name, PPCC = "ppccNormSinglyCensored", 
        `Shapiro-Wilk` = "swSinglyCensoredGofTestStatistic", 
        `Log-Likelihood` = "loglike.norm.boxcoxSinglyCensored")
    if (!optimize) {
        optimize.bounds <- rep(NA, 2)
        n <- length(lambda)
        objective.vec <- numeric(n)
        for (i in 1:n) {
            y <- boxcoxTransform(x = x, lambda = lambda[i], eps = eps)
            if (objective.name != "Log-Likelihood") {
                arg.list <- list(x = y, censored = censored, 
                  censoring.side = censoring.side)
            }
            else {
                est.list <- enormSinglyCensored(y, censored = censored, 
                  censoring.side = censoring.side, method = "mle")
                mean <- est.list$parameters["mean"]
                sd <- est.list$parameters["sd"]
                arg.list <- list(x = x, censored = censored, 
                  censoring.side = censoring.side, lambda = lambda[i], 
                  mean = mean, sd = sd, eps = eps)
            }
            objective.vec[i] <- do.call(objective.fcn, arg.list)
        }
    }
    else {
        fcn.to.min <- function(lambda, x, censored, censoring.side, 
            objective.name, objective.fcn, eps) {
            y <- boxcoxTransform(x = x, lambda = lambda, eps = eps)
            if (objective.name != "Log-Likelihood") {
                arg.list <- list(x = y, censored = censored, 
                  censoring.side = censoring.side)
            }
            else {
                est.list <- enormSinglyCensored(y, censored = censored, 
                  censoring.side = censoring.side, method = "mle")
                mean <- est.list$parameters["mean"]
                sd <- est.list$parameters["sd"]
                arg.list <- list(x = x, censored = censored, 
                  censoring.side = censoring.side, lambda = lambda, 
                  mean = mean, sd = sd, eps = eps)
            }
            -do.call(objective.fcn, arg.list)
        }
        optimize.bounds <- lambda
        nlminb.list <- nlminb(start = 1, objective = fcn.to.min, 
            lower = optimize.bounds[1], upper = optimize.bounds[2], 
            x = x, censored = censored, censoring.side = censoring.side, 
            objective.name = objective.name, objective.fcn = objective.fcn, 
            eps = eps)
        lambda <- nlminb.list$par
        objective.vec <- -nlminb.list$objective
    }
    names(optimize.bounds) <- c("lower", "upper")
    ret.list <- list(lambda = lambda, objective = objective.vec, 
        objective.name = objective.name, optimize = optimize, 
        optimize.bounds = optimize.bounds, eps = eps, sample.size = N, 
        censoring.side = censoring.side, censoring.levels = T1, 
        percent.censored = (100 * n.cen)/N, data.name = data.name, 
        censoring.name = censoring.name, bad.obs = bad.obs)
    if (include.x.and.censored) 
        ret.list <- append(ret.list, list(data = x, censored = censored), 
            after = 6)
    oldClass(ret.list) <- "boxcoxCensored"
    ret.list
}
