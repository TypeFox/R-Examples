`print.splm` <-
function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2,
                      quote = FALSE)
    } else {
        cat("No coefficients\n")
    }

    ## add printing of error variance parameters
    cat("\n")
    ec <- x$errcomp
    if (length(ec)) {
        cat("Error covariance parameters:\n")
        print.default(format(ec, digits = digits), print.gap = 2,
                      quote = FALSE)
    }

    else cat("No error covariance parameters\n")
    cat("\n")

    ## add printing of spatial autoregressive parameter
    ar <- x$arcoef
    if (length(ar)) {
        cat("\n")
        cat("Spatial autoregressive parameter:\n")
        print.default(format(ar, digits = digits), print.gap = 2,
                      quote = FALSE)
    }

    invisible(x)
}

