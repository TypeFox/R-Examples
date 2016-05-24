summary.pseudocpf <- function(object, conf.int = 0.95, scale = 1, ...) {
    if (!inherits(object, "pseudocpf")) {
        stop("'object' must be of class 'pseudocpf'")
    }
    nt <- length(object$timepoints)
    zzz <- list()
    zzz$call <- object$call
    beta <- object$fit$beta[-(1:nt)]
    se.beta <- sqrt(diag(object$fit$vbeta)[-(1:nt)])
    se.beta.ajs <- sqrt(diag(object$fit$vbeta.ajs)[-(1:nt)])
    if (all(se.beta.ajs == 0)) {
        tmp <- cbind(beta, exp(beta), se.beta, beta / se.beta,
                     1 - pchisq((beta/ se.beta)^2, 1))
    }
    else {
        tmp <- cbind(beta, exp(beta), se.beta.ajs, beta / se.beta.ajs,
                     1 - pchisq((beta/ se.beta.ajs)^2, 1))
    }
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", "se(coef)",
                                         "w", "Pr(>|w|)"))
    zzz$coefficients <- tmp
    z <- qnorm((1 + conf.int) / 2)
    tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se.beta),
                 exp(beta + z * se.beta))
    dimnames(tmp) <- list(names(beta),
                          c("exp(coef)", "exp(-coef)",
                            paste("lower .", round(100 * conf.int, 2), sep = ""),
                            paste("upper .", round(100 * conf.int, 2), sep = "")))
    zzz$conf.int <- tmp
    class(zzz) <- "summary.pseudocpf"
    zzz
}

print.summary.pseudocpf <- function(x, digits = max(getOption("digits") - 3, 3),
                                    signif.stars = getOption("show.signif.stars"),
                                    ...) {
    if (!inherits(x, "summary.pseudocpf")) {
        stop("'x' must be of class 'summary.pseudocpf'")
    }
    cat("Call:\n")
    dput(x$call)
    cat("\n")
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    printCoefmat(x$coefficients, digits = digits, signif.stars=signif.stars, ...)
    cat("\n")
    print(x$conf.int)
    cat("\n")
    invisible()
}
