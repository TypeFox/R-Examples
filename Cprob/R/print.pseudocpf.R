print.pseudocpf <- function(x, ...) {
    if (!inherits(x, "pseudocpf")) {
        stop("'x' must be of class 'pseudocpf'")
    }
    cat("Call:\n")
    dput(x$call); cat("\n")
    nt <- length(x$timepoints)
    coef <- x$fit$beta[-(1:nt)]
    se <- sqrt(diag(x$fit$vbeta)[-(1:nt)])
    se.ajs <- sqrt(diag(x$fit$vbeta.ajs))[-(1:nt)]
    if (all(se.ajs == 0)) {
        tmp <- cbind(coef, exp(coef), se, coef / se,
                     1 - pchisq((coef/ se)^2, 1))
    }
    else {
        tmp <- cbind(coef, exp(coef), se.ajs, coef / se.ajs,
                     1 - pchisq((coef/ se.ajs)^2, 1))
    }
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                         "se(coef)", "z", "p"))
    print(tmp)
    cat("\n")
    invisible()
}
