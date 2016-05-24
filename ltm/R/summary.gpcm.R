summary.gpcm <-
function (object, robust.se = FALSE, ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    coefs <- object$coefficients
    se <- sqrt(diag(vcov(object, robust = robust.se)))
    p <- length(coefs)
    ncatg <- sapply(coefs, length)
    constraint <- object$constraint
    spSE <- if (constraint == "gpcm") {
        ii <- rep(1:p, ncatg)
        split(se, ii)
    } else if (constraint == "1PL") {
        nse <- length(se)
        ii <- rep(1:p, ncatg - 1)
        lapply(split(se[-nse], ii), function (x) c(x, se[nse]))
    } else {
        ii <- rep(1:p, ncatg - 1)
        lapply(split(se, ii), function (x) c(x, NA))
    }
    coef.tab <- mapply(function (x, y) cbind(value = x, std.err = y, z.value = x / y), coefs, spSE, SIMPLIFY = FALSE)
    out <- list(coefficients = coef.tab)
    out$logLik <- object$log.Lik
    df <- sapply(object$coef, length)
    df <- switch(object$constraint,
        "gpcm" = sum(df),
        "1PL" = sum(df) - length(df) + 1,
        "rasch" = sum(df) - length(df)
    )
    out$AIC <- AIC(object)
    out$BIC <- AIC(object, k = log(attr(logLik(object), "n")))
    out$max.sc <- object$max.sc
    out$conv <- object$conv
    out$counts <- object$counts
    out$call <- object$call
    out$control <- object$control
    out$attr <- attr(object$X, "items")
    out$ancr <- attr(object$X, "anchoring")
    class(out) <- "summ.gpcm"
    out
}
