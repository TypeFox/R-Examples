`wa` <- function(x, ...) UseMethod("wa")

`wa.default` <-
    function(x, env,
             deshrink = c("inverse", "classical", "expanded", "none", "monotonic"),
             tol.dw = FALSE, useN2 = TRUE,
             na.tol = c("min","mean","max"),
             small.tol = c("min","mean","fraction","absolute"),
             min.tol = NULL, f = 0.1, ...)
{
    ## x = species abundances (weights), env = response vector
    x <- as.matrix(x)
    env <- as.numeric(env)
    ## drop species with no information
    if(any(csum <- colSums(x) == 0)) {
        x <- x[, !csum, drop = FALSE]
        warning("Some species contained no data. These have been deleted.")
    }
    ## drop samples with no species
    if(any(rsum <- rowSums(x) == 0)) {
        x <- x[!rsum, , drop = FALSE]
        env <- env[!rsum]
        warning("Some sites contained no data. These have been deleted.")
    }
    if(missing(deshrink))
        deshrink <- "inverse"
    deshrink <- match.arg(deshrink)
    if(missing(na.tol))
        na.tol <- "min"
    na.tol <- match.arg(na.tol)
    if(missing(small.tol))
        small.tol <- "min"
    small.tol <- match.arg(small.tol)
    ## do the WA
    fit <- waFit(x = x, y = env, tol.dw = tol.dw, useN2 = useN2,
                 deshrink = deshrink, na.tol = na.tol,
                 small.tol = small.tol, min.tol = min.tol, f = f)
    ## site/sample names need to be reapplied
    names(fit$fitted.values) <- rownames(x)
    ## species names need to be reapplied
    names(fit$wa.optima) <- colnames(x)
    ## residuals
    residuals <- env - fit$fitted.values
    ## RMSE of predicted/fitted values
    rmse <- sqrt(mean(residuals^2))
    ## r-squared
    r.squared <- cor(fit$fitted.values, env)^2
    ## bias statistics
    avg.bias <- mean(residuals)
    max.bias <- maxBias(residuals, env)
    ## the function call
    .call <- match.call()
    ## need to reset due to method dispatch
    .call[[1]] <- as.name("wa")
    ## returned object
    res <- list(wa.optima = fit$wa.optima,
                tolerances = fit$tolerances,
                model.tol = fit$model.tol,
                fitted.values = fit$fitted.values,
                residuals = residuals,
                coefficients = fit$coefficients,
                rmse = rmse, r.squared = r.squared,
                avg.bias = avg.bias, max.bias = max.bias,
                n.samp = fit$n.samp,
                n.spp = fit$n.spp,
                deshrink = deshrink, tol.dw = tol.dw,
                call = .call,
                orig.x = x, orig.env = env,
                options.tol =
                list(useN2 = useN2,
                     na.tol = na.tol,
                     small.tol = small.tol,
                     min.tol = min.tol,
                     f = f))
    class(res) <- "wa"
    res
}
