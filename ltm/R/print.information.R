print.information <-
function (x, ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Total Information =", round(x$InfoTotal, 2))
    cat("\nInformation in (", x$range[1], ", ", x$range[2], ") = ", round(x$InfoRange, 2), 
            " (", round(100 * x$InfoRange / x$InfoTotal, 2), "%)", sep = "")
    cat("\nBased on", if (is.null(x$items)) "all the items" else paste("items", paste(x$items, collapse = ", ")))
    cat("\n\n")
    invisible(x)
}
