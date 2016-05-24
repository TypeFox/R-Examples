## Print method for groc objects:
print.groc <- function(x, ...) {
    cat("Generalized Regression on Orthogonal Components optimized with the greed algorithm.")
    if (!is.null(x$validation))
        cat("\nCross-validated using", length(x$validation$segments),
            attr(x$validation$segments, "type"), "segments.")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    invisible(x)
}

## Summary method for groc objects
summary.groc <- function(object, what = "validation",
                        digits = 4, print.gap = 2, ...)
{
    if (is.null(object$validation)) what <- ""

    nobj <- nrow(object$scores)
    nresp <- length(object$Ymeans)
    yvarnames <- respnames(object)
    cat("Data: \tX dimension:", nobj, length(object$Xmeans),
        "\n\tY dimension:", nobj, nresp)
    cat("\nFit method:", object$method)
    cat("\nNumber of components considered:", object$ncomp,"\n")

        if (what == "validation") {
            cat("\nCross-validated using", length(object$validation$segments),
                attr(object$validation$segments, "type"), "segments.")
            cat("\nVALIDATION: PRESS\n")
            print(object$validation$PRESS, digits = digits, print.gap = print.gap, ...)
            cat("\nVALIDATION: PREMAD\n")
            print(object$validation$PREMAD, digits = digits, print.gap = print.gap, ...)
            cat("\nVALIDATION: RMSEP\n")
            print(object$validation$RMSEP, digits = digits, print.gap = print.gap, ...)
        }
    
}

