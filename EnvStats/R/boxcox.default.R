boxcox.default <-
function (x, lambda = {
    if (optimize) 
        c(-2, 2)
    else seq(-2, 2, by = 0.5)
}, optimize = FALSE, objective.name = "PPCC", eps = .Machine$double.eps, 
    include.x = TRUE, ...) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    if (any(x <= 0)) 
        stop("All non-missing, finite values of 'x' must be positive")
    if (!is.vector(lambda, mode = "numeric") || is.factor(lambda) || 
        any(!is.finite(lambda))) 
        stop("'lambda' must be a numeric vector with no missing or infinite values")
    if (optimize && length(unique(lambda)) != 2) 
        stop(paste("When optimize=TRUE, 'lambda' must be a vector", 
            "with two unique values that specify the lower and", 
            "upper bounds for the optimization"))
    if (optimize && (1 < min(lambda) || 1 > max(lambda))) 
        stop("When optimize=TRUE, the range of 'lambda' must contain 1")
    lambda <- sort(lambda)
    objective.name <- match.arg(objective.name, c("PPCC", "Shapiro-Wilk", 
        "Log-Likelihood"))
    objective.fcn <- switch(objective.name, PPCC = "ppccNorm", 
        `Shapiro-Wilk` = "swGofTestStatistic", `Log-Likelihood` = "loglike.norm.boxcox")
    sample.size <- length(x)
    if (!optimize) {
        optimize.bounds <- rep(NA, 2)
        n <- length(lambda)
        objective.vec <- numeric(n)
        for (i in 1:n) {
            y <- boxcoxTransform(x = x, lambda = lambda[i], eps = eps)
            if (objective.name != "Log-Likelihood") {
                arg.list <- list(x = y)
            }
            else {
                est.list <- enorm(y, method = "mle/mme")
                mean <- est.list$parameters["mean"]
                sd <- est.list$parameters["sd"]
                arg.list <- list(x = x, lambda = lambda[i], mean = mean, 
                  sd = sd, eps = eps)
            }
            objective.vec[i] <- do.call(objective.fcn, arg.list)
        }
    }
    else {
        fcn.to.min <- function(lambda, x, objective.name, objective.fcn, 
            eps) {
            y <- boxcoxTransform(x = x, lambda = lambda, eps = eps)
            if (objective.name != "Log-Likelihood") {
                arg.list <- list(x = y)
            }
            else {
                est.list <- enorm(y, method = "mle/mme")
                mean <- est.list$parameters["mean"]
                sd <- est.list$parameters["sd"]
                arg.list <- list(x = x, lambda = lambda, mean = mean, 
                  sd = sd, eps = eps)
            }
            -do.call(objective.fcn, arg.list)
        }
        optimize.bounds <- lambda
        nlminb.list <- nlminb(start = 1, objective = fcn.to.min, 
            lower = optimize.bounds[1], upper = optimize.bounds[2], 
            x = x, objective.name = objective.name, objective.fcn = objective.fcn, 
            eps = eps)
        lambda <- nlminb.list$par
        objective.vec <- -nlminb.list$objective
    }
    names(optimize.bounds) <- c("lower", "upper")
    ret.list <- list(lambda = lambda, objective = objective.vec, 
        objective.name = objective.name, optimize = optimize, 
        optimize.bounds = optimize.bounds, eps = eps, sample.size = sample.size, 
        data.name = data.name, bad.obs = bad.obs)
    if (include.x) 
        ret.list <- append(ret.list, list(data = x), after = 6)
    oldClass(ret.list) <- "boxcox"
    ret.list
}
