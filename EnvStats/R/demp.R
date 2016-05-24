demp <-
function (x, obs, discrete = FALSE, density.arg.list = NULL) 
{
    if (!is.vector(x, mode = "numeric") || (length.x <- length(x)) < 
        1) 
        stop("'x' must be a numeric vector with at least one element")
    names.x <- names(x)
    if (!is.vector(obs, mode = "numeric")) 
        stop("'obs' must be a numeric vector")
    if ((bad.obs <- sum(!(obs.ok <- is.finite(obs)))) > 0) {
        is.not.finite.warning(obs)
        obs <- obs[obs.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'obs' removed."))
    }
    if ((n.obs <- length(obs)) < 2) 
        stop("'obs' must be a numeric vector with at least 2 non-missing elements")
    if (!is.vector(discrete, mode = "logical") || length(discrete) != 
        1) 
        stop("'discrete' must be a logical scalar")
    na.index <- is.na(x)
    if (all(na.index)) 
        y <- as.numeric(rep(NA, length.x))
    else {
        y <- numeric(length.x)
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        x.no.na <- x[!na.index]
        if (discrete) {
            y.no.na <- table(obs)[match(x.no.na, sort(unique(obs)))]/n.obs
        }
        else {
            y.no.na <- approx(do.call("density", args = c(list(x = obs), 
                density.arg.list)), xout = x.no.na, method = "linear", 
                rule = 1)$y
        }
        y.no.na[is.na(y.no.na)] <- 0
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
