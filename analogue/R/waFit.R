`waFit` <- function(x, y, tol.dw, useN2, deshrink, na.tol,
                    small.tol, min.tol, f) {
    ## sample summaries
    n.samp <- nrow(x)
    n.spp <- ncol(x)
    ## calculate WA optima for each species in x
    wa.optima <- w.avg(x, y)
    ## fix-up NaN or NA optima
    if(any(miss <- is.na(wa.optima)))
        wa.optima[miss] <- 0
    ## compute tolerances
    tolerances <- tol <- w.tol(x, y, wa.optima, useN2 = useN2)
    ## fix-up tolerances for use in TF computations
    if(small.tol == "fraction") {
        if(!(f > 0 && f < 1))
            stop("'f' must be 0 < f < 1")
        frac <- f * diff(range(y))
        if(frac < min.tol)
            warning("Requested fraction of gradient is < minimum tolerance.")
    }
    tol <- fixUpTol(tol, na.tol = na.tol, small.tol = small.tol,
                    min.tol = min.tol, f = f, env = y)
    ## calculate WA estimate of env for each site
    wa.env <- if(tol.dw) {
        WATpred(x, wa.optima, tol, n.samp, n.spp)
    } else {
        WApred(x, wa.optima)
    }
    ## taken averages twice so deshrink
    expanded <- deshrink(y, wa.env, type = deshrink)
    wa.env <- expanded$env
    coefficients <- coef(expanded)
    ## returned object
    res <- list(wa.optima = wa.optima,
                tolerances = tolerances,
                model.tol = tol,
                fitted.values = wa.env,
                coefficients = coefficients,
                n.samp = n.samp,
                n.spp = n.spp)
    res
}
