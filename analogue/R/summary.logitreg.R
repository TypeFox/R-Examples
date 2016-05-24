`summary.logitreg` <- function(object, p = 0.9, ...) {
    FOO <- function(x, p) {
        coefs <- coef(summary(x))[2, ]
        dose <- dose.p(x, p = p)
        nobs <- length(x$y)
        IN <- sum(x$y)
        OUT <- sum(!x$y)
        c(IN, OUT, coefs, unname(dose[1]),
          unname(attr(dose, "SE")[,1]))
    }
    DF <- t(sapply(object$models, FOO, p = p, USE.NAMES = FALSE))
    DF <- data.frame(DF)
    names(DF) <- c("In","Out","E[Dij]","SE", "Z","p-value",
                   paste("Dij(p=", format(p), ")", sep = ""),
                   "SE (Dij)")
    class(DF) <- "summary.logitreg"
    DF
}

`print.summary.logitreg` <- function(x,
                                     digits = min(3, getOption("digits") - 4),
                                     ...) {
    class(x) <- "data.frame"
    x$`p-value` <- format.pval(x$`p-value`)
    want <- c(3:5, 7:8)
    x[, want] <- format(x[, want], digits = digits)
    cat("\n")
    writeLines(strwrap("Logit regression models"))
    cat("\n")
    print(x, ...)
    cat("\n")
    invisible(x)
}
