vcov.grm <-
function (object, ...) {
    if (!inherits(object, "grm"))
        stop("Use only with 'grm' objects.\n")
    if (is.null(object$hessian))
        stop("you should re-fit the model using 'Hessian = TRUE' in order to get standard errors.\n")
    inv.hes <- ginv(object$hessian)
    pars <- if (object$IRT.param) IRT.parm(object, digits.abbrv = object$control$digits.abbrv)$parms else object$coef
    p <- length(pars)
    ncatg <- sapply(pars, length)
    ind1 <- if (object$constrained) c(1, cumsum(ncatg[-p] - 1) + 1) else c(1, cumsum(ncatg[-p]) + 1)
    ind2 <- if (object$constrained) cumsum(ncatg - 1) else cumsum(ncatg)
    J <- jacobian(unlist(object$coef), object$constrained, ind1, ind2, p)
    for (i in seq(along = ind1)) {
        ind <- seq(ind1[i], if (object$constrained) ind2[i] else ind2[i] - 1)
        inv.hes[ind, ind] <- J[[i]] %*% inv.hes[ind, ind] %*% t(J[[i]])
    }
    nams <- if (object$constrained) {
        pars <- lapply(pars, function (x) x[-length(x)])
        c(paste(unlist(lapply(pars, names)), rep(names(pars), sapply(pars, length)), sep = "."), "Dscrmn")
    } else {
        paste(unlist(lapply(pars, names)), rep(names(pars), sapply(pars, length)), sep = ".")
    }
    dimnames(inv.hes) <- list(nams, nams)
    inv.hes
}
