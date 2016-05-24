print.GoF.rasch <-
function (x, ...) {
    if (!inherits(x, "GoF.rasch"))
        stop("Use only with 'GoF.rasch' objects.\n")
    p.val <- round(x$p.value, 3)
    p.val <- if (p.val < 0.001) "<0.001" else p.val
    cat("\nBootstrap Goodness-of-Fit using Pearson chi-squared\n")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Tobs:", round(x$Tobs, 2), "\n")
    cat("# data-sets:", x$B + 1, "\n")
    cat("p-value:", p.val, "\n\n")
    invisible(x)    
}
