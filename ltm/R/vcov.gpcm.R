vcov.gpcm <-
function (object, robust = FALSE, ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    inv.Hessian <- if (robust) {
        inv.H <- ginv(object$hessian)
        outer.score <- scoregpcmSNW(object)        
        inv.H %*% outer.score %*% inv.H
    } else
        ginv(object$hessian)
    nams <- if (object$constraint == "gpcm") {
        names(unlist(object$coefficients))
    } else if (object$constraint == "1PL") {
        nm <- lapply(object$coefficients, function (x) x[-length(x)])
        c(names(unlist(nm)), "alpha")
    } else {
        nm <- lapply(object$coefficients, function (x) x[-length(x)])
        names(unlist(nm))
    }
    dimnames(inv.Hessian) <- list(nams, nams)
    inv.Hessian
}
