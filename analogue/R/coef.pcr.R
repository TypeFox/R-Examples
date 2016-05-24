`coef.pcr` <- function(object, comps = NULL, ...) {
    coefs <- object$coefficients
    nc <- NCOL(coefs)
    if(is.null(comps))
        comps <- seq_len(nc)
    else {
        if(!is.numeric(comps))
            stop("Non-numeric selection of components requested.")
        if(min(comps) < 1 || max(comps) > nc)
            stop("Requested components outside range of actual components.")
    }
    coefs <- coefs[, comps]
    coefs
}
