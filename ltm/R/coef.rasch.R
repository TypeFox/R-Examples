coef.rasch <-
function (object, prob = FALSE, order = FALSE, ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    cof <- if (object$IRT.param){
            coefs <- IRT.parm(object)$parms
            p <- length(coefs)
            matrix(c(coefs[1:(p - 1)], rep(coefs[p], p - 1)), ncol = 2,
                    dimnames = list(colnames(object$X), c("Dffclt", "Dscrmn")))
        } else {
            object$coef
        }
    if (prob)
        cof <- cbind(cof, "P(x=1|z=0)" = plogis(object$coef[, 1]))
    if (order)
        cof <- cof[order(cof[, 1]), ]
    cof
}
