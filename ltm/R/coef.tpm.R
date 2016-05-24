coef.tpm <-
function (object, prob = FALSE, order = FALSE, ...) {
    if (!inherits(object, "tpm"))
        stop("Use only with 'tpm' objects.\n")
    cof <- if (object$IRT.param) IRT.parm(object)$parms else object$coef
    cof[, 1] <- plogis(cof[, 1]) * object$max.guessing
    if (prob)
        cof <- cbind(cof, "P(x=1|z=0)" = cof[, 1] + (1 - cof[, 1]) * plogis(object$coef[, 2]))
    if (order)
        cof <- cof[order(cof[, 2]), ]
    cof
}
