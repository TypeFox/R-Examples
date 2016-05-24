print.margins.ltm <-
function (x, digits = 2, ...) {
    if (!inherits(x, "margins.ltm"))
        stop("Use only with 'margins.ltm' objects.\n")    
    combs <- x$combs
    type <- x$type
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nFit on the", if (type == "two-way") "Two-Way" else "Three-Way", "Margins\n\n")
    for (i in 1:nrow(combs)) {
        cat("Response: (", paste(combs[i, ], collapse = ","), ")\n", sep = "")
        mat <- x$margins[, , i]
        mat <- mat[order(mat[, ncol(mat)], decreasing = TRUE), ]
        mat <- data.frame(round(mat[seq(1, x$nprint), ], digits))
        mat$rule <- ifelse(mat[, ncol(mat)] > x$rule, "***", "")
        names(mat) <- if (type == "two-way") {
            c("Item i", "Item j", "Obs", "Exp", "(O-E)^2/E", " ")
        } else {
            c("Item i", "Item j", "Item k", "Obs", "Exp", "(O-E)^2/E", " ")
        }
        print(mat)
        cat("\n")
    }
    if (any(apply(x$margins, 2, c)[, ncol(x$margins[, , 1])] > x$rule))
        cat("'***' denotes a chi-squared residual greater than", x$rule, "\n")
    cat("\n")
    invisible(x)
}
