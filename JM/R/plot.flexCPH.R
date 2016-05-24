plot.flexCPH <-
function (x, scale = c("survival", "cumHazard", "log-cumHazard"), survTimes = NULL, X = NULL, 
    xlab, lab, ylab, main, type, plot.it = TRUE, ...) {
    scale <- match.arg(scale)
    if (is.null(X))
        X <- x$X
    if (is.null(survTimes))
        survTimes <- seq(min(x$logT), max(x$logT), length.out = 51)
    W <- splineDesign(x$knots, survTimes, ord = x$control$ord)
    eta.w <- c(W %*% x$coefficients$gammas)
    eta.x <- c(X %*% x$coefficients$betas)
    eta <- outer(eta.x, eta.w, "+")
    fit <- switch(scale,
        "survival" = exp(- exp(eta)),
        "cumHazard" = exp(eta),
        "log-cumHazard" = eta)
    fit <- colMeans(fit)
    if (plot.it) {
        if (missing(xlab))
            xlab <- "Time"
        if (missing(ylab))
            ylab <- switch(scale,
                "survival" = "Survival",
                "cumHazard" = "Cumulative Hazard", 
                "log-cumHazard" = "log Cumulative Hazard")
        if (missing(main))
            main <- ""
        if (missing(type))
            type <- "l"
        plot(survTimes, fit, xlab = xlab, ylab = ylab, main = main, type = type, ...)
    }
    invisible(cbind(survTimes = exp(survTimes), fit = fit))
}
