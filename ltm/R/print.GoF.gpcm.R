print.GoF.gpcm <-
function (x, ...) {
    if (!inherits(x, "GoF.gpcm"))
        stop("Use only with 'GoF.gpcm' objects.\n")
    p.val <- round(x$p.value, 3)
    p.val <- if (p.val < 0.001) "<0.001" else p.val
     if (x$simulate.p.value) {
        cat("\nParametric Bootstrap Approximation to Pearson chi-squared Goodness-of-Fit Measure\n")
        cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
        cat("Tobs:", round(x$Tobs, 2), "\n")
        cat("# data-sets:", x$B + 1, "\n")
        cat("p-value:", p.val, "\n\n")
    } else {
        cat("\nPearson chi-squared Goodness-of-Fit Measure\n")
        cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
        cat("Tobs:", round(x$Tobs, 2), "\n")
        cat("df:", x$df, "\n")
        cat("p-value:", p.val, "\n\n")        
    }
    invisible(x)    
}
