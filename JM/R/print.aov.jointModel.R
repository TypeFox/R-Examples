print.aov.jointModel <-
function (x, ...) {
    if (!inherits(x, "aov.jointModel"))
        stop("Use only with 'aov.jointModel' objects.\n")
    if (is.null(x$L0)) {
        f <- function (dat) {
            dat[] <- lapply(dat, function (x) 
                round(unlist(x), 4))
            dat$'Pr(>|Chi|)' <- format.pval(dat$'Pr(>|Chi|)', 
                eps = 1e-04)
            dat
        }
        cat("\nMarginal Wald Tests Table\n")
        if (!is.null(x$aovTab.Y)) {
            cat("\nLongitudinal Process\n")
            print(f(x$aovTab.Y))
        }
        if (!is.null(x$aovTab.T)) {
            cat("\nEvent Process\n")
            print(f(x$aovTab.T))
        }
        if (!is.null(x$aovTab.L)) {
            cat("\nUser-defined Contrasts Matrix\n")
            print(f(x$aovTab.L))
        }
        cat("\n")
    } else {
        dat <- if (x$test) {
            p.val <- round(x$p.value, 4)
            p.val <- if (p.val < 0.0001) "<0.0001" else p.val
            data.frame(AIC = round(c(x$aic0, x$aic1), 2), 
                BIC = round(c(x$bic0, x$bic1), 2), 
                log.Lik = round(c(x$L0, x$L1), 2), 
                LRT = c(" ", round(x$LRT, 2)), df = c("", x$df), 
                p.value = c("", p.val), row.names = c(x$nam0, x$nam1))
        } else {
            data.frame(AIC = round(c(x$aic0, x$aic1), 2), 
                BIC = round(c(x$bic0, x$bic1), 2), 
                log.Lik = round(c(x$L0, x$L1), 2), df = c("", x$df), 
                row.names = c(x$nam0, x$nam1))    
        }
        cat("\n")
        print(dat)
        cat("\n")
    }
    invisible(x)
}
