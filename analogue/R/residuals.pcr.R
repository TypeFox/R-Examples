`residuals.pcr` <- function(object, comps = NULL, ...) {
    resi <- object$residuals
    nc <- NCOL(resi)
    if(is.null(comps))
        comps <- seq_len(nc)
    else {
        if(!is.numeric(comps))
            stop("Non-numeric selection of components requested.")
        if(min(comps) < 1 || max(comps) > nc)
            stop("Requested components outside range of actual components.")
    }
    resi <- resi[, comps]
    resi
}
