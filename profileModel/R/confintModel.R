`confintModel` <-
function (fitted, quantile = qchisq(0.95, 1), verbose = TRUE, 
    endpoint.tolerance = 0.001, max.zoom = 100, zero.bound = 1e-08, 
    stepsize = 0.5, stdn = 5, gridsize = 20, scale = FALSE, which = 1:length(coef(fitted)), 
    objective = stop("'objective' is missing."), agreement = TRUE, 
    method = "smooth", n.interpolations = 100, ...) 
{
    dotss <- match.call(expand.dots = FALSE)[["..."]]
    if (!(method %in% c("smooth", "zoom"))) 
        stop("Invalid method. The supported methods are 'smooth' and 'zoom'")
    prof <- profileModel(fitted, gridsize = gridsize, stdn = stdn, 
        quantile = quantile, objective = objective, agreement = agreement, 
        verbose = verbose, which = which, zero.bound = zero.bound, 
        stepsize = stepsize, scale = scale, ...)
    switch(method, zoom = ci <- profConfint(prof = prof, endpoint.tolerance = endpoint.tolerance, 
        max.zoom = max.zoom, verbose = verbose, method = "zoom"), 
        smooth = ci <- profConfint(prof = prof, n.interpolations = n.interpolations, 
            method = "smooth"))
    attr(ci, "profileModel object") <- NULL
    attr(ci, "fitted object") <- match.call()[["fitted"]]
    ci
}
