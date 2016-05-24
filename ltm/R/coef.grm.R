coef.grm <-
function (object, ...) {
    if (!inherits(object, "grm"))
        stop("Use only with 'grm' objects.\n")
    coefs <- if (object$IRT.param) IRT.parm(object, digits.abbrv = object$control$digits.abbrv)$parms else object$coef
    coefs <- lapply(coefs, round, digits = 3)
    if (all(sapply(coefs, length) == length(coefs[[1]])))
        coefs <- do.call(rbind, coefs)
    coefs
}
