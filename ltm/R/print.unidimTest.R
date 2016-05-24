print.unidimTest <-
function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nUnidimensionality Check using Modified Parallel Analysis\n")
    if (!is.null(x$call))
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nMatrix of tertachoric correlations\n")
    print(round(x$Rho, digits))
    cat("\nAlternative hypothesis: the second eigenvalue of the observed data is substantially larger",
        "\n\t\t\tthan the second eigenvalue of data under the assumed IRT model")
    cat("\n\nSecond eigenvalue in the observed data:", round(x$Tobs[2], digits))
    cat("\nAverage of second eigenvalues in Monte Carlo samples:", round(mean(x$T.boot[, 2], na.rm = TRUE), digits))
    cat("\nMonte Carlo samples:", NROW(x$T.boot))
    cat("\np-value:", round(x$p.value, digits))
    cat("\n\n")
    invisible(x)
}
