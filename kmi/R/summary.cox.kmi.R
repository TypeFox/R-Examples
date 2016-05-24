summary.cox.kmi <- function(object, conf.int = 0.95, scale = 1, ...) {
    if (!inherits(object, "cox.kmi")) {
        stop("'object' must be of class 'cox.kmi'")
    }
    if ("cox.penal" %in% class(object$cox.kmi.fit[[1]])) {
        stop("Doesn't work yet")
    }
    rval <- list()
    rval$call <- object$call
    beta <- object$coefficients
    se.beta <- sqrt(diag(object$var))
    tmp <- cbind(beta, exp(beta), se.beta, beta / se.beta,
                 2 * pt(abs(beta / se.beta), df = object$df, lower.tail = FALSE))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", "se(coef)",
                                         "t", "Pr(>|t|)"))
    rval$coefficients <- tmp
    z <- -qt((1 + conf.int)/2, object$df, lower.tail = FALSE)
    tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se.beta),
                 exp(beta + z * se.beta))
    dimnames(tmp) <- list(names(beta),
                          c("exp(coef)", "exp(-coef)",
                            paste("lower .", round(100 * conf.int, 2), sep = ""),
                            paste("upper .", round(100 * conf.int, 2), sep = "")))
    rval$conf.int <- tmp
    tmp <- lapply(object$individual.fit, summary, conf.int = conf.int, scale = scale, ...)
    rval$individual.fit <- tmp
    class(rval) <- "summary.cox.kmi"
    rval
}

