simulate.jointModel <-
function (object, nsim, seed = NULL, times = NULL, Data = NULL, ...) {
    thetas <- object$coefficients
    n <- object$n
    lag <- object$y$lag
    if (is.null(Data))
        Data <- object$data.id
    if (is.null(times))
        times <- sort(unique(object$times))
    method <- switch(object$method, "weibull-AFT-GH" = "weibull-AFT", 
        "weibull-PH-GH" = "weibull-PH", "spline-PH-GH" = "spline-PH",
        "piecewise-PH-GH" = "piecewise-PH", "Cox-PH-GH" =, 
        "ch-Laplace" = stop("not available.\n"))
    tt <- attr(delete.response(object$termsT), "term.labels")
    formulas <- list(Yfixed = reformulate(attr(delete.response(object$termsYx), 
        "term.labels")), Yrandom = object$formYz, timeVar = object$timeVar, 
        Tfixed = if (length(tt)) reformulate(tt) else reformulate("1"))
    simulateJM(nsim = nsim, nsub = n, thetas = thetas, times = times, 
        formulas = formulas, Data = Data, method = method, lag = lag, 
        seed = seed, ...)
}
