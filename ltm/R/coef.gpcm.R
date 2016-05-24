coef.gpcm <-
function (object, ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    coefs <- lapply(object$coefficients, round, digits = 3)
    if (all(sapply(coefs, length) == length(coefs[[1]])))
        coefs <- do.call(rbind, coefs)
    coefs
}
