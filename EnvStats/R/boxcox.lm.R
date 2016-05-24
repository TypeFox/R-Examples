boxcox.lm <-
function (x, lambda = {
    if (optimize) 
        c(-2, 2)
    else seq(-2, 2, by = 0.5)
}, optimize = FALSE, objective.name = "PPCC", eps = .Machine$double.eps, 
    include.x = TRUE, ...) 
{
    data.name <- deparse(substitute(x))
    if (is.null(x$y) || is.null(x$qr)) 
        x <- update(x, y = TRUE, qr = TRUE)
    Y <- x$y
    if (any(Y <= 0)) 
        stop("All non-missing values of the response variable must be positive")
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
        `Shapiro-Wilk` = "swGofTestStatistic", `Log-Likelihood` = stop("Log-Likelihood method not available"))
    sample.size <- length(Y)
    if (!optimize) {
        optimize.bounds <- rep(NA, 2)
        n <- length(lambda)
        objective.vec <- numeric(n)
        for (i in 1:n) {
            new.Y <- boxcoxTransform(x = Y, lambda = lambda[i], 
                eps = eps)
            residuals <- qr.resid(x$qr, new.Y)
            arg.list <- list(x = residuals)
            objective.vec[i] <- do.call(objective.fcn, arg.list)
        }
    }
    else {
        fcn.to.min <- function(lambda, x, Y, objective.fcn, eps) {
            new.Y <- boxcoxTransform(x = Y, lambda = lambda, 
                eps = eps)
            residuals <- qr.resid(x$qr, new.Y)
            arg.list <- list(x = residuals)
            -do.call(objective.fcn, arg.list)
        }
        optimize.bounds <- lambda
        nlminb.list <- nlminb(start = 1, objective = fcn.to.min, 
            lower = optimize.bounds[1], upper = optimize.bounds[2], 
            x = x, Y = Y, objective.fcn = objective.fcn, eps = eps)
        lambda <- nlminb.list$par
        objective.vec <- -nlminb.list$objective
    }
    names(optimize.bounds) <- c("lower", "upper")
    ret.list <- list(lambda = lambda, objective = objective.vec, 
        objective.name = objective.name, optimize = optimize, 
        optimize.bounds = optimize.bounds, eps = eps, sample.size = sample.size, 
        data.name = data.name)
    if (include.x) 
        ret.list <- append(ret.list, list(lm.obj = x), after = 6)
    oldClass(ret.list) <- "boxcoxLm"
    ret.list
}
